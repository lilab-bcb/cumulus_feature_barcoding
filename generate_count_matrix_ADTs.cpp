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
#include "ReadParser.hpp"

using namespace std;


const int totalseq_A_pos = 0;
const int totalseq_BC_pos = 10;


atomic<int> cnt, n_valid, n_valid_cell, n_valid_feature; // cnt: total number of reads; n_valid, reads with valid cell barcode and feature barcode; n_valid_cell, reads with valid cell barcode; n_valid_feature, reads with valid feature barcode

int n_threads, max_mismatch_cell, max_mismatch_feature, umi_len;
string feature_type, totalseq_type, scaffold_sequence;
int barcode_pos; // Antibody: Total-Seq A 0; Total-Seq B or C 10. Crispr: default 0, can be set by option
bool convert_cell_barcode;

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

	result_t(int cell_id, uint64_t umi, int feature_id) : cell_id(cell_id), feature_id(feature_id), umi(umi) {}
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
			inputs.push_back(move(one_pair));
		}

		input_dir = strtok(NULL, ",");
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


inline string safe_substr(const string& sequence, int pos, int length) {
	if (pos + length > sequence.length()) {
		printf("Error: Sequence length %d is too short (expected to be at least %d)!\n", (int)sequence.length(), pos + length);
		exit(-1);
	}
	return sequence.substr(pos, length);
}


inline bool extract_feature_barcode(const string& sequence, int feature_length, const string& feature_type, string& feature_barcode) {
	bool success = true;
	int start_pos, end_pos; // here start_pos and end_pos are with respect to feature sequence.

	if (feature_type == "antibody" || scaffold_sequence == "")
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


void detect_totalseq_type() {
	const int nskim = 10000; // Look at first 10000 reads.
	int ntotA, ntotBC, cnt;
	uint64_t binary_feature;
	Read read2;
	HashIterType feature_iter;

	cnt = ntotA = ntotBC = 0;
	for (auto&& input_pair : inputs) {
		iGZipFile gzip_in_r2(input_pair[1]);
		while (gzip_in_r2.next(read2) && cnt < nskim) {
			binary_feature = barcode_to_binary(safe_substr(read2.seq, totalseq_A_pos, feature_blen));
			feature_iter = feature_index.find(binary_feature);
			ntotA += (feature_iter != feature_index.end() && feature_iter->second.item_id >= 0);

			if (read2.seq.length() >= totalseq_BC_pos + feature_blen) {
				binary_feature = barcode_to_binary(safe_substr(read2.seq, totalseq_BC_pos, feature_blen));
				feature_iter = feature_index.find(binary_feature);
				ntotBC += (feature_iter != feature_index.end() && feature_iter->second.item_id >= 0);
			}
			++cnt;
		}
		if (cnt == nskim) break;
	}

	printf("ntotA = %d, ntotBC = %d.\n", ntotA, ntotBC);
	if (ntotA < 10 && ntotBC < 10) {
		printf("Error: Detected less than 10 feature barcodes in the first %d reads! Maybe you should consider to reverse complement your barcodes?\n", nskim);
		exit(-1);
	}

	totalseq_type = (ntotA > ntotBC ? "TotalSeq-A" : (umi_len == 12 ? "TotalSeq-B" : "TotalSeq-C"));
	barcode_pos = (totalseq_type == "TotalSeq-A" ? totalseq_A_pos : totalseq_BC_pos);
	printf("TotalSeq type is automatically detected as %s, barcode starts from 0-based position %d.\n", totalseq_type.c_str(), barcode_pos);
}


bool parse_feature_names(int n_feature, vector<string>& feature_names, int& n_cat, vector<string>& cat_names, vector<int>& cat_nfs, vector<int>& feature_categories) {
	std::size_t pos;
	string cat_str;

	pos = feature_names[0].find_first_of(',');
	if (pos == string::npos) {
		n_cat = 1;
		return false;
	}

	n_cat = 0;
	cat_names.clear();
	cat_nfs.clear();
	feature_categories.resize(n_feature, 0);
	for (int i = 0; i < n_feature; ++i) {
		pos = feature_names[i].find_first_of(',');
		assert(pos != string::npos);
		cat_str = feature_names[i].substr(pos + 1);
		feature_names[i] = feature_names[i].substr(0, pos);
		if (n_cat == 0 || cat_names.back() != cat_str) {
			cat_names.push_back(cat_str);
			cat_nfs.push_back(i);
			++n_cat;
		}
		feature_categories[i] = n_cat - 1;
	}
	cat_nfs.push_back(n_feature);

	return true;
}

void process_reads(ReadParser *parser, int thread_id) {
	string cell_barcode, umi, feature_barcode;
	uint64_t binary_cell, binary_umi, binary_feature;
	int read1_len;
	int cell_id, feature_id, collector_pos;
	bool valid_cell, valid_feature;
	HashIterType cell_iter, feature_iter;

	int cnt_, n_valid_, n_valid_cell_, n_valid_feature_;
	int delta_inc = 0; // counter for printing

	auto& buffer = result_buffer[thread_id];

	auto rg = parser->getReadGroup();
	while (parser->refill(rg)) {
		cnt_ = n_valid_ = n_valid_cell_ = n_valid_feature_ = 0;
		for (int i = 0; i < n_cat; ++i) buffer[i].clear();

		for (auto& read_pair : rg) {
			auto& read1 = read_pair[0];
			auto& read2 = read_pair[1];

			++cnt_;
			cell_barcode = safe_substr(read1.seq, 0, cell_blen);
			binary_cell = barcode_to_binary(cell_barcode);
			cell_iter = cell_index.find(binary_cell);
			valid_cell = cell_iter != cell_index.end() && cell_iter->second.item_id >= 0;

			valid_feature = extract_feature_barcode(read2.seq, feature_blen, feature_type, feature_barcode);
			if (valid_feature) {
				binary_feature = barcode_to_binary(feature_barcode);
				feature_iter = feature_index.find(binary_feature);
				valid_feature = feature_iter != feature_index.end() && feature_iter->second.item_id >= 0;				
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
				binary_umi = barcode_to_binary(umi);

				cell_id = cell_iter->second.item_id;
				feature_id = feature_iter->second.item_id;
				collector_pos = detected_ftype ? feature_categories[feature_id] : 0;
				buffer[collector_pos].emplace_back(cell_id, binary_umi, feature_id);
			}
		}

		for (int i = 0; i < n_cat; ++i) {
			auto& dataCollector = dataCollectors[i];
			collector_locks[i]->lock();
			for (auto& r : buffer[i]) dataCollector.insert(r.cell_id, r.umi, r.feature_id);
			collector_locks[i]->unlock();
		}

		cnt += cnt_;
		n_valid += n_valid_;
		n_valid_cell += n_valid_cell_;
		n_valid_feature += n_valid_feature_;

		delta_inc += cnt_;
		if (delta_inc >= 1000000) {
			printf("Processed %d reads.\n", cnt.load());
			delta_inc = 0;
		}
	}
}


int main(int argc, char* argv[]) {
	if (argc < 5) {
		printf("Usage: generate_count_matrix_ADTs cell_barcodes.txt[.gz] feature_barcodes.csv fastq_folders output_name [-p #] [--max-mismatch-cell #] [--feature feature_type] [--max-mismatch-feature #] [--umi-length len] [--barcode-pos #] [--convert-cell-barcode] [--scaffold-sequence sequence] [--chunk chunk_size]\n");
		printf("Arguments:\n\tcell_barcodes.txt[.gz]\t10x genomics barcode white list, either gzipped or not.\n");
		printf("\tfeature_barcodes.csv\tfeature barcode file;barcode,feature_name[,feature_category]. Optional feature_category is required only if hashing and citeseq data share the same sample index.\n");
		printf("\tfastq_folders\tfolder containing all R1 and R2 FASTQ files ending with 001.fastq.gz .\n");
		printf("\toutput_name\toutput file name prefix.\n");
		printf("Options:\n");
		printf("\t-p #\tnumber of threads. This number should be >= 2. [default: 2]\n");
		printf("\t--max-mismatch-cell #\tmaximum number of mismatches allowed for cell barcodes. [default: 0]\n");
		printf("\t--feature feature_type\tfeature type can be either antibody or crispr. [default: antibody]\n");
		printf("\t--max-mismatch-feature #\tmaximum number of mismatches allowed for feature barcodes. [default: 2]\n");
		printf("\t--umi-length len\tlength of the UMI sequence. [default: 12]\n");
		printf("\t--barcode-pos #\tstart position of barcode in read 2, 0-based coordinate. [default: automatically determined for antibody; 0 for crispr]\n");
		printf("\t--convert-cell-barcode\tconvert cell barcode to match RNA cell barcodes for 10x Genomics' data. Note that both cmo and 10x crispr need to set this option to convert feature barcoding barcodes to RNA barcodes. When data is hashing/CITE-Seq, this option will be automatically turned on for TotalSeq-B antibodies.\n");
		printf("\t--scaffold-sequence sequence\tscaffold sequence used to locate the protospacer for sgRNA. This option is only used for crispr data. If --barcode-pos is not set and this option is set, try to locate barcode in front of the specified scaffold sequence.\n");
		printf("Outputs:\n\toutput_name.csv\tfeature-cell count matrix. First row: [Antibody/CRISPR],barcode_1,...,barcode_n;Other rows: feature_name,feature_count_1,...,feature_count_n.\n");
		printf("\toutput_name.stat.csv.gz\tSufficient statistics file. First row: Barcode,UMI,Feature,Count; Other rows: each row describe the read count for one barcode-umi-feature combination.\n\n");
		printf("\tIf feature_category presents, this program will output the above two files for each feature_category. For example, if feature_category is hashing, output_name.hashing.csv and output_name.hashing.stat.csv.gz will be generated.\n");
		printf("\toutput_name.report.txt\tA report file summarizing barcode, UMI and read results.\n");
		exit(-1);
	}

	start_ = time(NULL);

	n_threads = 2;
	max_mismatch_cell = 1;
	feature_type = "antibody";
	max_mismatch_feature = 3;
	umi_len = 10;
	barcode_pos = -1;
	totalseq_type = "";
	scaffold_sequence = "";
	convert_cell_barcode = false;

	size_t chunk_size = 100000;

	for (int i = 5; i < argc; ++i) {
		if (!strcmp(argv[i], "-p")) {
			n_threads = atoi(argv[i + 1]);
		}
		if (!strcmp(argv[i], "--max-mismatch-cell")) {
			max_mismatch_cell = atoi(argv[i + 1]);
		}
		if (!strcmp(argv[i], "--feature")) {
			feature_type = argv[i + 1];
		}
		if (!strcmp(argv[i], "--max-mismatch-feature")) {
			max_mismatch_feature = atoi(argv[i + 1]);
		}
		if (!strcmp(argv[i], "--umi-length")) {
			umi_len = atoi(argv[i + 1]);
		}
		if (!strcmp(argv[i], "--barcode-pos")) {
			barcode_pos = atoi(argv[i + 1]);
		}
		if (!strcmp(argv[i], "--convert-cell-barcode")) {
			convert_cell_barcode = true;
		}
		if (!strcmp(argv[i], "--scaffold-sequence")) {
			scaffold_sequence = argv[i + 1];
		}

		if (!strcmp(argv[i], "--chunk")) {
			chunk_size = atoi(argv[i+1]);
		}

	}

	printf("Load feature barcodes.\n");
	parse_sample_sheet(argv[2], n_feature, feature_blen, feature_index, feature_names, max_mismatch_feature);
	detected_ftype = parse_feature_names(n_feature, feature_names, n_cat, cat_names, cat_nfs, feature_categories);

	parse_input_directory(argv[3]);

	if (feature_type == "antibody") {
		if (barcode_pos < 0) detect_totalseq_type(); // if specify --barcode-pos, must be a customized assay
	}
    else {
		if (feature_type != "crispr") {
			printf("Do not support unknown feature type %s!\n", feature_type.c_str());
			exit(-1);
		}
		if (barcode_pos < 0) barcode_pos = 0; // default is 0
	}

	interim_ = time(NULL);
	printf("Load cell barcodes.\n");
	convert_cell_barcode = convert_cell_barcode || (feature_type == "antibody" && totalseq_type == "TotalSeq-B");
	parse_sample_sheet(argv[1], n_cell, cell_blen, cell_index, cell_names, max_mismatch_cell, convert_cell_barcode);
	end_ = time(NULL);
	printf("Time spent on parsing cell barcodes = %.2fs.\n", difftime(end_, interim_));

	interim_ = end_;

	int np, nt;

	np = 1 ? n_threads <= 4 || inputs.size() == 1 : 2;
	nt = n_threads - np;

	dataCollectors.resize(n_cat);
	result_buffer.resize(nt);
	for (int i = 0; i < nt; ++i) result_buffer[i].resize(n_cat);
	for (int i = 0; i < n_cat; ++i) collector_locks.emplace_back(new mutex());

	cnt = 0;
	n_valid = 0;
	n_valid_cell =0 ;
	n_valid_feature = 0;

	ReadParser *parser = new ReadParser(inputs, nt, np, chunk_size);

	for (int i = 0; i < nt; ++i)
		processingThreads_.emplace_back([parser, i](){ process_reads(parser, i); });

	for (auto& thread : processingThreads_) thread.join();
	delete parser;
	result_buffer.clear();

	end_ = time(NULL);
	printf("Parsing input data is finished. Time spent = %.2fs.\n", difftime(end_, interim_));
	interim_ = end_;

	string output_name = argv[4];
	ofstream fout;

	fout.open(output_name + ".report.txt");
	fout<< "Total number of reads: "<< cnt<< endl;
	fout<< "Number of reads with valid cell barcodes: "<< n_valid_cell<< " ("<< fixed<< setprecision(2)<< n_valid_cell * 100.0 / cnt << "%)"<< endl;
	fout<< "Number of reads with valid feature barcodes: "<< n_valid_feature<< " ("<< fixed<< setprecision(2)<< n_valid_feature * 100.0 / cnt << "%)"<< endl;
	fout<< "Number of reads with valid cell and feature barcodes: "<< n_valid<< " ("<< fixed<< setprecision(2)<< n_valid * 100.0 / cnt << "%)"<< endl;

	if (n_cat == 0)
		dataCollectors[0].output(output_name, feature_type, 0, n_feature, cell_names, umi_len, feature_names, fout, n_threads);
	else
		for (int i = 0; i < n_cat; ++i) {
			printf("Feature '%s':\n", cat_names[i].c_str());
			dataCollectors[i].output(output_name + "." + cat_names[i], feature_type, cat_nfs[i], cat_nfs[i + 1], cell_names, umi_len, feature_names, fout, n_threads);
		}
	fout.close();

	printf("%s.report.txt is written.\n", output_name.c_str());
	
	end_ = time(NULL);
	printf("Outputs are written. Time spent = %.2fs", difftime(end_, interim_));
	
	printf("Total time spent (not including destruct objects) = %.2fs.\n", difftime(end_, start_));

	return 0;
}
