#include <ctime>
#include <cstdio>
#include <cstdint>
#include <cassert>
#include <cstring>
#include <cstdlib>
#include <string>
#include <vector>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <memory>
#include <atomic>
#include <mutex>
#include <thread>

#include "dirent.h"

#include "gzip_utils.hpp"
#include "barcode_utils.hpp"
#include "datamatrix_utils.hpp"
#include "auto_utils.hpp"
#include "ReadParser.hpp"

using namespace std;


string cb_dir;
string chemistry;

atomic<int> cnt, n_valid, n_valid_cell, n_valid_feature, n_reads_valid_umi, prev_cnt; // cnt: total number of reads; n_valid, reads with valid cell barcode and feature barcode; n_valid_cell, reads with valid cell barcode; n_valid_feature, reads with valid feature barcode; prev_cnt: for printing # of reads processed purpose

int n_threads, max_mismatch_cell, max_mismatch_feature, umi_len;
float read_ratio_cutoff;
int umi_count_cutoff;
string genome, feature_type, totalseq_type, scaffold_sequence, umi_correct_method;
int barcode_pos; // Antibody: Total-Seq A 0; Total-Seq B or C 10. Crispr: default 0, can be set by option

time_t start_, interim_, end_;

vector<vector<string>> inputs;

int n_cell, n_feature; // number of cell and feature barcodes
int cell_blen, feature_blen; // cell barcode length and feature barcode length
vector<string> cell_names, feature_names;
HashType cell_index, feature_index;

int n_cat; // number of feature categories (e.g. hashing, citeseq)
bool detected_ftype; // if feature csv contains feature type information
vector<string> cat_names; // category names
vector<int> cat_nfs, feature_categories; // cat_nfs, number of features in each category; int representing categories.
vector<DataCollector> dataCollectors;

struct result_t {
	int cell_id, feature_id;
	uint64_t umi;

	result_t(int cell_id, int feature_id, uint64_t umi) : cell_id(cell_id), feature_id(feature_id), umi(umi) {}
};

vector<vector<vector<result_t>>> result_buffer;
vector<thread> processingThreads_;
vector<unique_ptr<mutex>> collector_locks;


void parse_input_directory(char* input_dirs) {
	DIR *dir;
	struct dirent *ent;
	vector<string> mate1s, mate2s;

	string mate1_pattern = string("R1_001.fastq.gz");
	string mate2_pattern = string("R2_001.fastq.gz");
	string dir_name;

	char *input_dir = strtok(input_dirs, ",");

	inputs.clear();
	while (input_dir != NULL) {
		assert((dir = opendir(input_dir)) != NULL);

		dir_name = string(input_dir) + "/";

		mate1s.clear();
		mate2s.clear();

		while ((ent = readdir(dir)) != NULL) {
			if (ent->d_type == DT_REG) {
				string file_name = string(ent->d_name);
				size_t pos;

				pos = file_name.find(mate1_pattern);
				if (pos != string::npos && pos + mate1_pattern.length() == file_name.length()) {
					mate1s.push_back(file_name);
				}

				pos = file_name.find(mate2_pattern);
				if (pos != string::npos && pos + mate2_pattern.length() == file_name.length()) {
					mate2s.push_back(file_name);
				}
			}
		}

		int s = mate1s.size();

		assert(s == mate2s.size());
		sort(mate1s.begin(), mate1s.end());
		sort(mate2s.begin(), mate2s.end());

		for (int i = 0; i < s; ++i) {
			vector<string> one_pair(2);
			one_pair[0] = dir_name + mate1s[i];
			one_pair[1] = dir_name + mate2s[i];
			inputs.push_back(std::move(one_pair));
		}

		input_dir = strtok(NULL, ",");
	}

	if (inputs.empty()) {
		printf("No FASTQ file found in input folder(s): \"%s\"!\n", input_dirs);
		exit(-1);
	}
}


// return rightmost position + 1
inline int matching(const string& readseq, const string& pattern, int nmax_mis, int pos, int& best_value) {
	int nmax_size = nmax_mis * 2 + 1;
	int f[2][7]; // for banded dynamic programming, max allowed mismatch = 3
	// f[x][y] : x, pattern, y, readseq
	// f[x][y] = min(f[x - 1][y - 1] + delta, f[x][y - 1] + 1, f[x - 1][y] + 1)
	int rlen = readseq.length(), plen = pattern.length();
	int prev, curr, rpos;
	int value, best_j = -1;

	// init f[-1], do not allow insertion at the beginning
	for (int j = 0; j < nmax_size; ++j) f[1][j] = nmax_mis + 1;
	f[1][nmax_mis] = 0;

	// Dynamic Programming
	prev = 1; curr = 0;
	best_value = 0;
	int i;
	for (i = 0; i < plen; ++i) {
		best_value = nmax_mis + 1; best_j = -1;
		for (int j = 0; j < nmax_size; ++j) {
			value = nmax_mis + 1;
			rpos = pos + i + (j - nmax_mis);
			if (rpos >= 0 && rpos < rlen) value = min(value, f[prev][j] + (pattern[i] != readseq[rpos])); // match/mismatch
			if (j > 0) value = min(value, f[curr][j - 1] + 1); // insertion
			if (j + 1 < nmax_size) value = min(value, f[prev][j + 1] + 1); // deletion
			f[curr][j] = value;
			if (best_value > value) { best_value = value; best_j = j; }
		}
		if (best_value > nmax_mis) break;
		prev = curr; curr ^= 1;
	}
	return best_value <= nmax_mis ? pos + i + (best_j - nmax_mis) : -1;
}


// [start, end]
inline int locate_scaffold_sequence(const string& sequence, const string& scaffold, int start, int end, int max_mismatch) {
	int i, pos, best_value = max_mismatch + 1, value;

	for (i = start; i <= end; ++i) {
		pos = matching(sequence, scaffold, max_mismatch, i, best_value);
		if (pos >= 0) break;
	}

	if (best_value > 0) {
		for (int j = i + 1; j <= i + max_mismatch; ++j) {
			pos = matching(sequence, scaffold, max_mismatch, j, value);
			if (best_value > value) best_value = value, i = j;
		}
	}

	return i <= end ? i : -1;
}


inline bool extract_feature_barcode(const string& sequence, int feature_length, const string& feature_type, string& feature_barcode) {
	bool success = true;
	int start_pos, end_pos; // here start_pos and end_pos are with respect to feature sequence.

	if (feature_type != "crispr" || scaffold_sequence == "")
		feature_barcode = safe_substr(sequence, barcode_pos, feature_length);
	else {
		// With scaffold sequence, locate it first
		start_pos = 0;
		end_pos = locate_scaffold_sequence(sequence, scaffold_sequence, start_pos + feature_length - max_mismatch_feature, sequence.length() - (scaffold_sequence.length() - 2), 2);
		success = end_pos >= 0;
		if (success) {
			if (end_pos - start_pos >= feature_length)
				feature_barcode = safe_substr(sequence, end_pos - feature_length, feature_length);
			else
				feature_barcode = string(feature_length - (end_pos - start_pos), 'N') + safe_substr(sequence, start_pos, end_pos - start_pos);
		}
	}
	return success;
}

void process_reads(ReadParser *parser, int thread_id) {
	string cell_barcode, umi, feature_barcode;
	uint64_t binary_cell, binary_umi, binary_feature;
	int read1_len;
	int cell_id, feature_id, collector_pos;
	bool valid_cell, valid_feature;
	HashIterType cell_iter, feature_iter;

	int cnt_, n_valid_, n_valid_cell_, n_valid_feature_, n_reads_valid_umi_;

	auto& buffer = result_buffer[thread_id];

	auto rg = parser->getReadGroup();
	while (parser->refill(rg)) {
		cnt_ = n_valid_ = n_valid_cell_ = n_valid_feature_ = n_reads_valid_umi_ = 0;
		for (int i = 0; i < n_cat; ++i) buffer[i].clear();

		for (auto& read_pair : rg) {
			auto& read1 = read_pair[0];
			auto& read2 = read_pair[1];

			++cnt_;
			cell_barcode = safe_substr(read1.seq, 0, cell_blen);
			binary_cell = barcode_to_binary(cell_barcode);
			cell_iter = cell_index.find(binary_cell);
			valid_cell = cell_iter != cell_index.end() && cell_iter->second.vid >= 0;

			valid_feature = extract_feature_barcode(read2.seq, feature_blen, feature_type, feature_barcode);
			if (valid_feature) {
				binary_feature = barcode_to_binary(feature_barcode);
				feature_iter = feature_index.find(binary_feature);
				valid_feature = feature_iter != feature_index.end() && feature_iter->second.vid >= 0;
			}

			n_valid_cell_ += valid_cell;
			n_valid_feature_ += valid_feature;

			if (valid_cell && valid_feature) {
				++n_valid_;
				read1_len = read1.seq.length();
				if (read1_len < cell_blen + umi_len) {
					printf("Warning: Processing thread %d detected read1 length %d is smaller than cell barcode length %d + UMI length %d. Shorten UMI length to %d!\n", thread_id, read1_len, cell_blen, umi_len, read1_len - cell_blen);
					umi_len = read1_len - cell_blen;
				}
				umi = safe_substr(read1.seq, cell_blen, umi_len);
				binary_umi = barcode_to_binary(umi, true);
				if (binary_umi == INVALID_UMI) continue;
				++n_reads_valid_umi_;

				cell_id = cell_iter->second.vid;
				feature_id = feature_iter->second.vid;
				collector_pos = detected_ftype ? feature_categories[feature_id] : 0;
				buffer[collector_pos].emplace_back(cell_id, feature_id, binary_umi);
			}
		}

		for (int i = 0; i < n_cat; ++i) {
			auto& dataCollector = dataCollectors[i];
			collector_locks[i]->lock();
			for (auto& r : buffer[i]) dataCollector.insert(r.cell_id, r.feature_id, r.umi);
			collector_locks[i]->unlock();
		}

		cnt += cnt_;
		n_valid += n_valid_;
		n_valid_cell += n_valid_cell_;
		n_valid_feature += n_valid_feature_;
		n_reads_valid_umi += n_reads_valid_umi_;

		if (cnt - prev_cnt >= 1000000) {
			printf("Processed %d reads.\n", cnt.load());
			prev_cnt = cnt.load();
		}
	}
}


int main(int argc, char* argv[]) {
	if (argc < 5) {
		printf("Usage: generate_feature_barcode_matrices cell_barcodes_dir feature_type feature_barcodes.csv fastq_folders output_name [OPTIONS]\n");
		printf("Arguments:\n\tcell_barcodes_dir\tPath to the folder containing 10x genomics barcode inclusion list files.\n");
		printf("\tfeature_type\tfeature type can be: \'hashing\', \'citeseq\', \'cmo\', \'crispr\', \'adt\' (both hashing and citeseq features in the same sample).\n");
		printf("\tfeature_barcodes.csv\tfeature barcode file;barcode,feature_name[,feature_category]. Optional feature_category is required only if hashing and citeseq data share the same sample index.\n");
		printf("\tfastq_folders\tfolder containing all R1 and R2 FASTQ files ending with \'001.fastq.gz\'.\n");
		printf("\toutput_name\toutput file name prefix.\n");
		printf("Options:\n");
		printf("\t-p #\tnumber of threads. This number should be >= 2. [default: 3]\n");
		printf("\t--genome genome_name\tgenome reference name. [default: \'\']\n");
		printf("\t--chemistry chemistry_type\tchemistry type. [default: auto]\n");
		printf("\t--max-mismatch-feature #\tmaximum number of mismatches allowed for feature barcodes. [default: 2]\n");
		printf("\t--barcode-pos #\tstart position of barcode in read 2, 0-based coordinate. [default: automatically determined for antibody; 0 for crispr]\n");
		printf("\t--scaffold-sequence sequence\tscaffold sequence used to locate the protospacer for sgRNA. This option is only used for crispr data. If --barcode-pos is not set and this option is set, try to locate barcode in front of the specified scaffold sequence.\n");
		printf("\t--umi-correct-method method_name\tUMI correction method to use. Available options: \'cluster\', \'adjacency\', \'directional\'. [default: directional]\n");
		printf("\t--umi-count-cutoff #\tRead count threshold (non-inclusive) to filter UMIs. Only works when <feature_type> is \'crispr\'.  [default: 0]\n");
		printf("\t--read-ratio-cutoff #\tRead count ratio threshold (non-inclusive) to filter chimeric reads. Only works when <feature_type> is \'crispr\'  [default: 0.5]\n");
		printf("Outputs:\n\t<output_name>.<feature_type>.h5\tfeature-cell count matrix in 10x HDF5 format.\n");
		printf("\t<output_name>.<feature_type>.molecule_info.h5\tSufficient statistics file in a simplied 10x HDF5 format. Each entry describe the read count for one barcode-umi-feature combination.\n\n");
		printf("\tIf feature_category presents, this program will output the above two files for each feature_category. For example, if feature_category is hashing, <output_name>.hashing.h5 and <output_name>.hashing.molecule_info.h5 will be generated.\n");
		printf("\t<output_name>.<feature_type>.report.txt\tA report file summarizing barcode, UMI and read results.\n");
		exit(-1);
	}

	start_ = time(NULL);

	n_threads = 3;
	genome = "";
	chemistry = "auto";
	max_mismatch_cell = -1;
	max_mismatch_feature = 2;
	umi_len = -1;
	umi_correct_method = "directional";
	umi_count_cutoff = 0;
	read_ratio_cutoff = 0.5;
	barcode_pos = -1;
	totalseq_type = "";
	scaffold_sequence = "";

	for (int i = 5; i < argc; ++i) {
		if (!strcmp(argv[i], "-p")) {
			n_threads = stoi(argv[i + 1]);
		}
		if (!strcmp(argv[i], "--genome")) {
			genome = argv[i + 1];
		}
		if (!strcmp(argv[i], "--chemistry")) {
			chemistry = argv[i + 1];
		}
		if (!strcmp(argv[i], "--max-mismatch-feature")) {
			max_mismatch_feature = stoi(argv[i + 1]);
		}
		if (!strcmp(argv[i], "--umi-correct-method")) {
			umi_correct_method = argv[i + 1];
		}
		if (!strcmp(argv[i], "--umi-count-cutoff")) {
			umi_count_cutoff = stoi(argv[i + 1]);
		}
		if (!strcmp(argv[i], "--read-ratio-cutoff")) {
			read_ratio_cutoff = stof(argv[i + 1]);
			if (read_ratio_cutoff >= 1.0 || read_ratio_cutoff < 0) {
				printf("\'--read-ratio-cutoff\' must be a fractional number within [0, 1)!\n");
				exit(-1);
			}
		}
		if (!strcmp(argv[i], "--barcode-pos")) {
			barcode_pos = stoi(argv[i + 1]);
		}
		if (!strcmp(argv[i], "--scaffold-sequence")) {
			scaffold_sequence = argv[i + 1];
		}
	}

	feature_type = argv[2];

	printf("Load feature barcodes.\n");
	parse_sample_sheet(argv[3], n_feature, feature_blen, feature_index, feature_names, max_mismatch_feature);
	detected_ftype = parse_feature_names(n_feature, feature_index, feature_names, n_cat, cat_names, cat_nfs, feature_categories);

	parse_input_directory(argv[4]);

	cb_dir = argv[1];
	if (cb_dir.length() > 0 && cb_dir[cb_dir.length()-1] != '/')
		cb_dir += "/";

	// Determine chemistry, totalseq_type (for antibody assays), barcode_pos, umi_len, max_mismatch_cells
	auto_detection(cb_dir, inputs, feature_type, feature_blen, feature_index, umi_len, chemistry, totalseq_type, max_mismatch_cell, barcode_pos, scaffold_sequence);

	interim_ = time(NULL);
	printf("Load cell barcodes.\n");
	parse_sample_sheet(cb_inclusion_file_dict[chemistry], n_cell, cell_blen, cell_index, cell_names, max_mismatch_cell);
	end_ = time(NULL);
	printf("Time spent on parsing cell barcodes = %.2fs.\n", difftime(end_, interim_));

	interim_ = end_;

	int np = min(max(1, n_threads / 3), (int)inputs.size());
	int nt = np * 2;

	dataCollectors.resize(n_cat);
	result_buffer.resize(nt);
	for (int i = 0; i < nt; ++i) result_buffer[i].resize(n_cat);
	for (int i = 0; i < n_cat; ++i) collector_locks.emplace_back(new mutex());

	cnt = 0; prev_cnt = 0;
	n_valid = 0;
	n_valid_cell =0 ;
	n_valid_feature = 0;
	n_reads_valid_umi = 0;

	ReadParser *parser = new ReadParser(inputs, nt, np);

	for (int i = 0; i < nt; ++i)
		processingThreads_.emplace_back([parser, i](){ process_reads(parser, i); });

	for (auto& thread : processingThreads_) thread.join();
	delete parser;
	result_buffer.clear();

	end_ = time(NULL);
	printf("Parsing input data is finished. %d reads are processed. Time spent = %.2fs.\n", cnt.load(), difftime(end_, interim_));
	interim_ = end_;

	string output_name = argv[5];
	ofstream fout;

	fout.open(output_name + "." + feature_type + ".report.txt");
	fout<< "Total number of reads: "<< cnt<< endl;
	fout<< "Number of reads with valid cell barcodes: "<< n_valid_cell<< " ("<< fixed<< setprecision(2)<< n_valid_cell * 100.0 / cnt << "%)"<< endl;
	fout<< "Number of reads with valid feature barcodes: "<< n_valid_feature<< " ("<< fixed<< setprecision(2)<< n_valid_feature * 100.0 / cnt << "%)"<< endl;
	fout<< "Number of reads with valid cell and feature barcodes: "<< n_valid<< " ("<< fixed<< setprecision(2)<< n_valid * 100.0 / cnt << "%)"<< endl;
	fout<< "Number of reads with valid cell, feature and UMI barcodes: "<< n_reads_valid_umi<< " ("<< fixed<< setprecision(2)<< n_reads_valid_umi * 100.0 / cnt << "%)" << endl;

	int total_umis_raw, total_cells_raw, cur_total_umis, cur_total_cells;
	string cur_ftype;
	int feature_start, feature_end;
	double runtime_umi_correct, runtime_write_output;
	runtime_umi_correct = runtime_write_output = 0;
	time_t start_1, end_1;

	for (int i = 0; i < n_cat; ++i) {
		if (!detected_ftype) {
			cur_ftype = feature_type;
			feature_start = 0;
			feature_end = n_feature;
		} else {
			cur_ftype = cat_names[i];
			feature_start = cat_nfs[i];
			feature_end = cat_nfs[i + 1];
		}

		// Generate raw count matrix
		start_1 = time(NULL);
		dataCollectors[i].output(output_name + "." + cur_ftype, "raw", genome, cur_ftype, feature_start, feature_end, cell_names, umi_len, feature_names, fout, n_threads, false);
		end_1 = time(NULL);
		runtime_write_output += difftime(end_1, start_1);
		total_umis_raw = dataCollectors[i].get_total_umis();
		total_cells_raw = dataCollectors[i].get_total_cells();

		// UMI correction
		if (i == 0)
			printf("[UMI correction] Use %s method for correction.\n", umi_correct_method.c_str());
		start_1 = time(NULL);
		dataCollectors[i].correct_umi(umi_len, umi_correct_method);
		end_1 = time(NULL);
		runtime_umi_correct += difftime(end_1, start_1);
		cur_total_umis = dataCollectors[i].get_total_umis();
		printf("[UMI correction] After UMI correction, %d (%.2f%%) UMIs are kept.\n", cur_total_umis, cur_total_umis * 1.0 / total_umis_raw * 100);

		// Generate UMI correct count matrix
		bool verbose_report = cur_ftype != "crispr" ? true : false;
		start_1 = time(NULL);
		dataCollectors[i].output(output_name + "." + cur_ftype, "umi_correct", genome, cur_ftype, feature_start, feature_end, cell_names, umi_len, feature_names, fout, n_threads, verbose_report);
		end_1 = time(NULL);
		runtime_write_output += difftime(end_1, start_1);

		// Chimeric filtering if needed
		if (cur_ftype == "crispr") {
			if (umi_count_cutoff > 0)  printf("[Chimeric filter] UMI count filtering by cutoff %d.\n", umi_count_cutoff);
			else  printf("[Chimeric filter] No UMI count filtering.\n");

			if (std::abs(read_ratio_cutoff) < std::numeric_limits<float>::epsilon())
				printf("[Chimeric filter] No PCR chimeric filtering.\n");
			else {
				printf("[Chimeric filter] PCR chimeric filtering by ratio cutoff %.2f.\n", read_ratio_cutoff);
				start_1 = time(NULL);
				dataCollectors[i].filter_chimeric_reads(umi_count_cutoff, read_ratio_cutoff);
				end_1 = time(NULL);
				runtime_umi_correct += difftime(end_1, start_1);
				cur_total_cells = dataCollectors[i].get_total_cells();
				cur_total_umis = dataCollectors[i].get_total_umis();
				printf("[Chimeric filter] After PCR chimeric filtering, %d (%.2f%%) cells and %d (%.2f%%) UMIs are kept.\n",
					cur_total_cells, cur_total_cells * 1.0 / total_cells_raw * 100,
					cur_total_umis, cur_total_umis * 1.0 / total_umis_raw * 100
				);
			}

			// Generate chimeric filter count matrix
			start_1 = time(NULL);
			dataCollectors[i].output(output_name + "." + cur_ftype, "chimeric_filter", genome, cur_ftype, feature_start, feature_end, cell_names, umi_len, feature_names, fout, n_threads, true);
			end_1 = time(NULL);
			runtime_write_output += difftime(end_1, start_1);
		}
	}
	printf("UMI correction is finished. Time spent = %.2fs.\n", runtime_umi_correct);
	printf("Outputs are written. Time spent = %.2fs.\n", runtime_write_output);

	fout.close();
	printf("%s.%s.report.txt is written.\n", output_name.c_str(), feature_type.c_str());
	end_ = time(NULL);

	printf("Total time spent (not including destruct objects) = %.2fs.\n", difftime(end_, start_));

	return 0;
}
