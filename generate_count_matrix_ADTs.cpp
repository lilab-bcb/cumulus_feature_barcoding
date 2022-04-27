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

#include "dirent.h"

#include "gzip_utils.hpp"
#include "barcode_utils.hpp"
#include "datamatrix_utils.hpp"

using namespace std;


const int STRLEN = 1005;

const int totalseq_A_pos = 0;
const int totalseq_BC_pos = 10;

struct InputFile{
	string input_r1, input_r2;

	InputFile(string r1, string r2) : input_r1(r1), input_r2(r2) {}
};

int max_mismatch_cell, max_mismatch_feature, umi_len;
string feature_type, totalseq_type, scaffold_sequence;
int barcode_pos; // Antibody: Total-Seq A 0; Total-Seq B or C 10. Crispr: default 0, can be set by option
bool convert_cell_barcode;

time_t start_time, end_time;

vector<InputFile> inputs; 

Read read1, read2;

int n_cell, n_feature; // number of cell and feature barcodes
int cell_blen, feature_blen; // cell barcode length and feature barcode length
vector<string> cell_names, feature_names;
HashType cell_index, feature_index;
HashIterType cell_iter, feature_iter;

int f[2][7]; // for banded dynamic programming, max allowed mismatch = 3


int n_cat; // number of feature categories (e.g. hashing, citeseq)
vector<string> cat_names; // category names
vector<int> cat_nfs, feature_categories; // cat_nfs, number of features in each category; int representing categories.
vector<DataCollector> dataCollectors;



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
			inputs.emplace_back(dir_name + mate1s[i], dir_name + mate2s[i]);
		}

		input_dir = strtok(NULL, ",");
	}
}


// return rightmost position + 1
inline int matching(const string& readseq, const string& pattern, int nmax_mis, int pos, int& best_value) {
	int nmax_size = nmax_mis * 2 + 1;
	// f[x][y] : x, pattern, y, readseq
	// f[x][y] = min(f[x - 1][y - 1] + delta, f[x][y - 1] + 1, f[x - 1][y] + 1)
	int rlen = readseq.length(), plen = pattern.length();
	int prev, curr, rpos;
	int value, best_j;

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
	int i, pos, best_value, value;

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
	int start_pos, end_pos, best_value; // here start_pos and end_pos are with respect to feature sequence.

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

	cnt = ntotA = ntotBC = 0;
	for (auto&& input_fastq : inputs) {
		iGZipFile gzip_in_r2(input_fastq.input_r2);
		while (gzip_in_r2.next(read2) == 4 && cnt < nskim) {
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


void parse_feature_names(int n_feature, vector<string>& feature_names, int& n_cat, vector<string>& cat_names, vector<int>& cat_nfs, vector<int>& feature_categories) {
	std::size_t pos;
	string cat_str;

	n_cat = 0;

	pos = feature_names[0].find_first_of(',');
	if (pos != string::npos) {
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
	}
}


int main(int argc, char* argv[]) {
	if (argc < 5) {
		printf("Usage: generate_count_matrix_ADTs cell_barcodes.txt feature_barcodes.csv fastq_folders output_name [--max-mismatch-cell #] [--feature feature_type] [--max-mismatch-feature #] [--umi-length len] [--barcode-pos #] [--convert-cell-barcode] [--scaffold-sequence sequence]\n");
		printf("Arguments:\n\tcell_barcodes.txt\t10x genomics barcode white list\n");
		printf("\tfeature_barcodes.csv\tfeature barcode file;barcode,feature_name[,feature_category]. Optional feature_category is required only if hashing and citeseq data share the same sample index\n");
		printf("\tfastq_folders\tfolder contain all R1 and R2 FASTQ files ending with 001.fastq.gz\n");
		printf("\toutput_name\toutput file name prefix\n");
		printf("Options:\n\t--max-mismatch-cell #\tmaximum number of mismatches allowed for cell barcodes [default: 1]\n");
		printf("\t--feature feature_type\tfeature type can be either antibody or crispr [default: antibody]\n");
		printf("\t--max-mismatch-feature #\tmaximum number of mismatches allowed for feature barcodes [default: 3]\n");
		printf("\t--umi-length len\tlength of the UMI sequence [default: 10]\n");
		printf("\t--barcode-pos #\tstart position of barcode in read 2, 0-based coordinate [default: automatically determined for antibody; 0 for crispr].\n");
		printf("\t--convert-cell-barcode\tconvert cell barcode to match RNA cell barcodes for 10x Genomics' data. Note that both cmo and 10x crispr need to set this option to convert feature barcoding barcodes to RNA barcodes. When data is hashing/CITE-Seq, this option will be automatically turned on for TotalSeq-B antibodies.\n");
		printf("\t--scaffold-sequence sequence\tscaffold sequence used to locate the protospacer for sgRNA. This option is only used for crispr data. If --barcode-pos is not set and this option is set, try to locate barcode in front of the specified scaffold sequence.\n");
		printf("Outputs:\n\toutput_name.csv\tfeature-cell count matrix. First row: [Antibody/CRISPR],barcode_1,...,barcode_n;Other rows: feature_name,feature_count_1,...,feature_count_n\n");
		printf("\toutput_name.stat.csv\tSufficient statistics file. First row: Barcode,UMI,Feature,Count; Other rows: each row describe the read count for one barcode-umi-feature combination\n\n");
		printf("\tIf feature_category presents, this program will output the above two files for each feature_category. For example, if feature_category is hashing, output_name.hashing.csv and output_name.hashing.stat.csv.gz will be generated.\n");
		printf("\toutput_name.report.txt\tA report file summarizing barcode, UMI and read results.\n");
		exit(-1);
	}

	start_time = time(NULL);

	max_mismatch_cell = 1;
	feature_type = "antibody";
	max_mismatch_feature = 3;
	umi_len = 10;
	barcode_pos = -1;
	totalseq_type = "";
	scaffold_sequence = "";
	convert_cell_barcode = false;

	for (int i = 5; i < argc; ++i) {
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
	}

	printf("Load feature barcodes.\n");
	parse_sample_sheet(argv[2], n_feature, feature_blen, feature_index, feature_names, max_mismatch_feature);
	parse_feature_names(n_feature, feature_names, n_cat, cat_names, cat_nfs, feature_categories);

	parse_input_directory(argv[3]);

	if (feature_type == "antibody") {
		if (barcode_pos < 0) detect_totalseq_type(); // if specify --barcode-pos, must be a customized assay
	} else {
		if (feature_type != "crispr") {
			printf("Do not support unknown feature type %s!\n", feature_type.c_str());
			exit(-1);
		}
		if (barcode_pos < 0) barcode_pos = 0; // default is 0
	}

	printf("Load cell barcodes.\n");
	convert_cell_barcode = convert_cell_barcode || (feature_type == "antibody" && totalseq_type == "TotalSeq-B");
	parse_sample_sheet(argv[1], n_cell, cell_blen, cell_index, cell_names, max_mismatch_cell, convert_cell_barcode);
	printf("Time spent on parsing cell barcodes = %.2fs.\n", difftime(time(NULL), start_time));

	int cnt = 0;
	string cell_barcode, umi, feature_barcode;
	uint64_t binary_cell, binary_umi, binary_feature;
	int read1_len;
	int feature_id, collector_pos;

	dataCollectors.resize(n_cat > 0 ? n_cat : 1);

	for (auto&& input_fastq : inputs) {
		iGZipFile gzip_in_r1(input_fastq.input_r1.c_str());
		iGZipFile gzip_in_r2(input_fastq.input_r2.c_str());
		while (gzip_in_r1.next(read1) == 4 && gzip_in_r2.next(read2) == 4) {
			++cnt;
			
			cell_barcode = safe_substr(read1.seq, 0, cell_blen);
			binary_cell = barcode_to_binary(cell_barcode);
			cell_iter = cell_index.find(binary_cell);

			if (cell_iter != cell_index.end() && cell_iter->second.item_id >= 0) {
				if (extract_feature_barcode(read2.seq, feature_blen, feature_type, feature_barcode)) {
					binary_feature = barcode_to_binary(feature_barcode);
					feature_iter = feature_index.find(binary_feature);
					if (feature_iter != feature_index.end() && feature_iter->second.item_id >= 0) {
						read1_len = read1.seq.length();
						if (read1_len < cell_blen + umi_len) {
							printf("Warning: Detected read1 length %d is smaller than cell barcode length %d + UMI length %d. Shorten UMI length to %d!\n", read1_len, cell_blen, umi_len, read1_len - cell_blen);
							umi_len = read1_len - cell_blen;
						}
						umi = safe_substr(read1.seq, cell_blen, umi_len);
						binary_umi = barcode_to_binary(umi);

						feature_id = feature_iter->second.item_id;
						collector_pos = n_cat > 0 ? feature_categories[feature_id] : 0;
						dataCollectors[collector_pos].insert(cell_iter->second.item_id, binary_umi, feature_id);
					}
				}
			}

			if (cnt % 1000000 == 0) printf("Processed %d reads.\n", cnt);
		}
	}

	printf("Parsing input data is finished.\n");

	string output_name = argv[4];

	int n_valid = 0;
	ofstream fout;
	fout.open(output_name + ".report.txt");
	fout<< "Total number of reads: "<< cnt<< endl<< endl;

	if (n_cat == 0)
		n_valid = dataCollectors[0].output(output_name, feature_type, 0, n_feature, cell_names, umi_len, feature_names, fout);
	else 
		for (int i = 0; i < n_cat; ++i) {
			printf("Feature '%s':\n", cat_names[i].c_str());
			n_valid += dataCollectors[i].output(output_name + "." + cat_names[i], feature_type, cat_nfs[i], cat_nfs[i + 1], cell_names, umi_len, feature_names, fout);
		}

	fout<< "Number of reads with valid cell and feature barcodes: "<< n_valid<< " ("<< fixed<< setprecision(2)<< n_valid * 100.0 / cnt << "%)"<< endl;
	fout.close();

	end_time = time(NULL);
	printf("Time spent = %.2fs.\n", difftime(end_time, start_time));

	return 0;
}
