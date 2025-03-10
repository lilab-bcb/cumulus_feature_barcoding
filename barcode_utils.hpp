#ifndef BARCODE_UTILS
#define BARCODE_UTILS

#include <cstdio>
#include <cassert>
#include <cstdint>
#include <string>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <numeric>

struct ValueType {
	int item_id;
	char n_mis; // number of mismatches

	ValueType() : item_id(-1), n_mis(0) {}
	ValueType(int item_id, char n_mis) : item_id(item_id), n_mis(n_mis) {}
};

typedef std::unordered_map<uint64_t, ValueType> HashType;
typedef HashType::iterator HashIterType;

const int STEP = 3;
const int BASE = 7;
const int UPPER = 21;
const int NNUC = 5; // ACGTN

const char id2base[NNUC] = {'A', 'C', 'G', 'T', 'N'};

// Assume char's range is -128..127
const int CHAR_RANGE = 128;

static std::vector<int> init_base2id() {
	std::vector<int> vec(CHAR_RANGE, -1);
	vec['a'] = vec['A'] = 0;
	vec['c'] = vec['C'] = 1;
	vec['g'] = vec['G'] = 2;
	vec['t'] = vec['T'] = 3;
	vec['n'] = vec['N'] = 4;

	return vec;
}

static const std::vector<int> base2id = init_base2id();

static std::vector<char> init_base2rcbase() {
	std::vector<char> vec(CHAR_RANGE, -1);
	vec['A'] = 'T'; vec['C'] = 'G'; vec['G'] = 'C'; vec['T'] = 'A'; vec['N'] = 'N';
	return vec;
}

static const std::vector<char> base2rcbase = init_base2rcbase();

static std::vector<std::vector<uint64_t> > init_aux_arr() {
	std::vector<std::vector<uint64_t> > aux_arr;
	for (int i = 0; i < UPPER; ++i) {
		std::vector<uint64_t> arr(NNUC + 1, 0);
		for (uint64_t j = 0; j < NNUC; ++j) arr[j] = j << (STEP * i);
		arr[NNUC] = uint64_t(BASE) << (STEP * i);
		aux_arr.push_back(arr);
	}
	return aux_arr;
}

static const std::vector<std::vector<uint64_t> > aux_arr = init_aux_arr();

uint64_t barcode_to_binary(const std::string& barcode) {
	uint64_t binary_id = 0;
	char c;
	if (barcode.length() > UPPER) {
		printf("Barcode %s exceeds the length limit %d!\n", barcode.c_str(), UPPER);
		exit(-1);
	}
	for (auto&& it = barcode.rbegin(); it != barcode.rend(); ++it) {
		c = *it;
		if (base2id[c] < 0) {
			printf("Barcode %s contains unknown bases %c!\n", barcode.c_str(), c);
			exit(-1);
		}
		binary_id <<= STEP;
		binary_id += base2id[c];
	}
	return binary_id;
}

std::string binary_to_barcode(uint64_t binary_id, int len) {
	std::string barcode(len, 0);
	for (int i = 0; i < len; ++i) {
		barcode[i] = id2base[binary_id & BASE];
		binary_id >>= STEP;
	}
	return barcode;
}

inline bool insert(HashType& index_dict, uint64_t key, ValueType&& value) {
	std::pair<HashIterType, bool> ret;
	ret = index_dict.insert(std::make_pair(key, value));
	if (ret.second) return true;
	if (ret.first->second.n_mis == 0 && value.n_mis == 0) {
		printf("Cumulus identified two identical barcodes! Please check your barcode file.\n");
		exit(-1);
	}
	if (ret.first->second.n_mis == 0 || value.n_mis == 0) {
		printf("Mismatch value is too large. Please decrease the number of mismatches allowed.\n");
		exit(-1);
	}
	ret.first->second.item_id = -1;
	return false;
}

inline void mutate_index_one_mismatch(HashType& index_dict, std::string& barcode, int item_id) {
	int len = barcode.length();
	uint64_t binary_id = barcode_to_binary(barcode);

	insert(index_dict, binary_id, ValueType(item_id, 0));
	for (int i = 0; i < len; ++i) {
		uint64_t val = binary_id & aux_arr[i][NNUC];
		for (int j = 0; j < NNUC; ++j)
			if (val != aux_arr[i][j]) {
				insert(index_dict, binary_id - val + aux_arr[i][j], ValueType(item_id, 1));
			}
	}
}

inline void mutate_index(HashType& index_dict, uint64_t binary_id, int len, int item_id, int max_mismatch, int mismatch, int pos) {
	insert(index_dict, binary_id, ValueType(item_id, mismatch));
	if (mismatch >= max_mismatch) return;

	for (int i = pos; i < len; ++i) {
		uint64_t val = binary_id & aux_arr[i][NNUC];
		for (int j = 0; j < NNUC; ++j)
			if (val != aux_arr[i][j]) {
				mutate_index(index_dict, binary_id - val + aux_arr[i][j], len, item_id, max_mismatch, mismatch + 1, i + 1);
			}
	}
}

inline void ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {
        return !std::isspace(ch);
    }));
}

inline void rtrim(std::string& s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
        return !std::isspace(ch);
    }).base(), s.end());
}

// trim from both ends (in place)
inline void trim(std::string &s) {
    rtrim(s);
    ltrim(s);
}

inline void group_by_modality(HashType& index_dict, std::vector<std::string>& index_names) {
	std::vector<int> indices(index_names.size());
	std::iota(indices.begin(), indices.end(), 0);
	std::sort(indices.begin(), indices.end(),
		[&index_names](int l, int r) {
			std::string s1 = index_names[l];
			std::string s2 = index_names[r];
			return s1.substr(s1.find_first_of(',') + 1) < s2.substr(s2.find_first_of(',') + 1);
		}
	);

	bool already_sorted = true;
	for (int i = 0; i < indices.size(); ++i)
		if (indices[i] != i) {
			already_sorted = false;
			break;
		}

	// No action if barcodes are already grouped by modality column
	if (already_sorted) return;

	std::vector<int> idx_map(indices.size(), -1);
	std::vector<std::string> tmp_names(index_names);
	for (int i = 0; i < indices.size(); ++i) {
		idx_map[indices[i]] = i;
		index_names[i] = tmp_names[indices[i]];
	}
	for (auto iter = index_dict.begin(); iter != index_dict.end(); ++iter) {
		iter->second.item_id = idx_map[iter->second.item_id];
	}
}

inline void parse_one_line(const std::string& line, int& n_barcodes, int& barcode_len, HashType& index_dict, std::vector<std::string>& index_names, int max_mismatch) {
	std::string index_name, index_seq;
	std::size_t pos;

	if (line.empty()) return;

	std::string delimiters = ",\t";    // Feature barcode file uses ','; Tow-column cell barcode file uses '\t'.
	pos = line.find_first_of(delimiters);

	if (pos != std::string::npos) { index_seq = line.substr(0, pos); trim(index_seq); index_name = line.substr(pos + 1); trim(index_name); }    // Feature barcode file, or two-column cell barcode file
	else { index_seq = line; trim(index_seq); index_name = index_seq; }    // One-column cell barcode file

	if (index_seq.empty() && index_name.empty()) return;

	if (barcode_len == 0) barcode_len = index_seq.length();
	else assert(barcode_len == index_seq.length());

	if (max_mismatch == 1) mutate_index_one_mismatch(index_dict, index_seq, n_barcodes);
	else mutate_index(index_dict, barcode_to_binary(index_seq), index_seq.length(), n_barcodes, max_mismatch, 0, 0);

	index_names.push_back(index_name);
	++n_barcodes;
}

inline void skip_bom(std::string& line) {
	size_t start = 0;

	if (line.length() >= 3 && line.substr(0, 3) == "\xEF\xBB\xBF")   // UTF-8
		start = 3;
	else if (line.length() >= 2 && (line.substr(0, 2) == "\xFF\xFE" or line.substr(0, 2) == "\xFE\xFF"))    // UTF-16
		start = 2;
	else if (line.length() >= 4 && (line.substr(0, 4) == "\x00\x00\xFE\xFF" or line.substr(0, 4) == "\xFF\xFE\x00\x00"))    // UTF-32
		start = 4;

	line = line.substr(start);
}

void parse_sample_sheet(const std::string& sample_sheet_file, int& n_barcodes, int& barcode_len, HashType& index_dict, std::vector<std::string>& index_names, int max_mismatch = 1, bool verbose = true) {
	std::string line;

	n_barcodes = 0;
	barcode_len = 0;
	index_dict.clear();
	index_names.clear();

	bool is_first_line = true;

	if (sample_sheet_file.length() > 3 && sample_sheet_file.substr(sample_sheet_file.length() - 3, 3) == ".gz") { // input sample sheet is gzipped
		iGZipFile gin(sample_sheet_file);
		while (gin.next(line)) {
			if (is_first_line) {
				skip_bom(line);
				is_first_line = false;
			}
			parse_one_line(line, n_barcodes, barcode_len, index_dict, index_names, max_mismatch);
		}
	}
	else {
		std::ifstream fin(sample_sheet_file);
		while (std::getline(fin, line)) {
			if (is_first_line) {
				skip_bom(line);
				is_first_line = false;
			}
			parse_one_line(line, n_barcodes, barcode_len, index_dict, index_names, max_mismatch);
		}
		fin.close();
	}
	if (verbose) printf("%s is parsed. n_barcodes = %d, and barcode_len = %d.\n", sample_sheet_file.c_str(), n_barcodes, barcode_len);

	int n_amb = 0;
	for (auto&& kv : index_dict)
		if (kv.second.item_id < 0) ++n_amb;
	if (verbose) printf("In the index, %d out of %d items are ambigious, percentage = %.2f%%.\n", n_amb, (int)index_dict.size(), n_amb * 100.0 / index_dict.size());
}

#endif
