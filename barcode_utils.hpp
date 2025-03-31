#ifndef BARCODE_UTILS
#define BARCODE_UTILS

#include <cstdio>
#include <cassert>
#include <cstdint>
#include <string>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <queue>
#include <algorithm>
#include <numeric>
#include <bit>
#include <utility>


struct BinaryCodeType {
	uint64_t bid; // binary ID
	uint32_t mask; // mask is 1 if the corresponding position is N or mismatch

	BinaryCodeType(): bid(0), mask(0) {}
	BinaryCodeType(uint64_t bid, uint32_t mask): bid(bid), mask(mask) {}
};

struct ValueType {
	int vid;
	uint32_t mask; // mask indicates location of mismatches

	ValueType() : vid(-1), mask(0) {}
	ValueType(int vid, uint32_t mask) : vid(vid), mask(mask) {}

	bool is_mutated_from(uint32_t& m) {
		return (mask & m) == m;
	}
};

struct IndexType {
	uint64_t bid;
	int vid;
	uint32_t mask;

	IndexType(): bid(0), vid(0), mask(0) {}
	IndexType(uint64_t bid, int vid, uint32_t mask) : bid(bid), vid(vid), mask(mask) {}
};

typedef std::unordered_map<uint64_t, ValueType> HashType;
typedef HashType::iterator HashIterType;

const int STEP = 2;
const int BASE = 3;
const int UPPER = 32;
const int NNUC = 4; // ACGT

const char id2base[NNUC] = {'A', 'C', 'G', 'T'};

// Assume char's range is -128..127
const int CHAR_RANGE = 128;

static std::vector<int> init_base2id() {
	std::vector<int> vec(CHAR_RANGE, -1);
	vec['a'] = vec['A'] = 0;
	vec['c'] = vec['C'] = 1;
	vec['g'] = vec['G'] = 2;
	vec['t'] = vec['T'] = 3;
	return vec;
}

static const std::vector<int> base2id = init_base2id();

static std::vector<char> init_base2rcbase() {
	std::vector<char> vec(CHAR_RANGE, -1);
	vec['A'] = 'T'; vec['C'] = 'G'; vec['G'] = 'C'; vec['T'] = 'A';
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

BinaryCodeType barcode_to_binary(const std::string& barcode) {
	char c;
	BinaryCodeType binary_code;

	if (barcode.length() > UPPER) {
		printf("Barcode %s exceeds the length limit %d!\n", barcode.c_str(), UPPER);
		exit(-1);
	}

	for (auto&& it = barcode.rbegin(); it != barcode.rend(); ++it) {
		c = *it;
		if (base2id[c] < 0 && c != 'N') {
			printf("Barcode %s contains unknown bases %c!\n", barcode.c_str(), c);
			exit(-1);
		}
		binary_code.bid <<= STEP;
		binary_code.mask <<= 1;
		if (base2id[c] >= 0)
			binary_code.bid += base2id[c];
		else
			binary_code.mask += 1;
	}

	return binary_code;
}

std::string binary_to_barcode(BinaryCodeType binary_code, int len) {
	std::string barcode(len, 0);
	for (int i = 0; i < len; ++i) {
		if (binary_code.mask & 1)
			barcode[i] = 'N';
		else
			barcode[i] = id2base[binary_code.bid & BASE];
		binary_code.mask >>= 1;
		binary_code.bid >>= STEP;
	}
	return barcode;
}

inline bool insert(HashType& index_dict, uint64_t key, ValueType&& value, IndexType& track) {
	std::pair<HashIterType, bool> ret;
	ret = index_dict.insert(std::make_pair(key, value));
	track.bid = key;
	track.vid = value.vid;
	track.mask = value.mask;
	if (ret.second) return true;
	if (ret.first->second.mask == 0 && value.mask == 0) {
		printf("Cumulus identified two identical barcodes! Please check your barcode file.\n");
		exit(-1);
	}
	if (std::popcount(ret.first->second.mask) == std::popcount(value.mask)) {
		ret.first->second.vid = -1;
		return false;
	}
	return true;
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
		iter->second.vid = idx_map[iter->second.vid];
	}
}

inline void parse_one_line(const std::string& line, int& n_barcodes, int& barcode_len, HashType& index_dict, std::vector<std::string>& index_names, int max_mismatch, std::queue<IndexType>* buffer) {
	std::string index_name, index_seq;
	std::size_t pos;

	if (line.empty()) return;

	std::string delimiters = ",\t";    // Feature barcode file uses ','; Tow-column cell barcode file uses '\t'.
	pos = line.find_first_of(delimiters);

	if (pos != std::string::npos) { index_seq = line.substr(0, pos); trim(index_seq); index_name = line.substr(pos + 1); trim(index_name); }    // Feature barcode file, or two-column cell barcode file
	else { index_seq = line; trim(index_seq); index_name = index_seq; }    // One-column cell barcode file

	if (index_seq.empty() && index_name.empty()) return;

	if (barcode_len == 0) {
		barcode_len = index_seq.length();
		assert(barcode_len <= UPPER);
	}
	else
		assert(barcode_len == index_seq.length());

	BinaryCodeType bin_code = barcode_to_binary(index_seq);
	IndexType index_record;
	insert(index_dict, bin_code.bid, ValueType(n_barcodes, 0), index_record);    // assume N not in cell/feature barcodes
	buffer->emplace(index_record);

	index_names.emplace_back(index_name);
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

void insert_mutated_index(HashType& index_dict, int barcode_len, int max_mismatch, std::queue<IndexType>* buffer1, std::queue<IndexType>* buffer2, bool verbose = true) {
	int cur_mismatch = 1;
	while (cur_mismatch <= max_mismatch) {
		while (!buffer1->empty()) {
			IndexType& index_val = buffer1->front();
			if (std::popcount(index_val.mask) < max_mismatch) {
				// Start from the next position since the last mutation
				int start_pos = UPPER - std::countl_zero(index_val.mask);
				for (int i = start_pos; i < barcode_len; ++i) {
					uint64_t val = index_val.bid & aux_arr[i][NNUC];
					// Set the position in mask to 1
					uint32_t mask = index_val.mask | (1 << i);
					for (int j = 0; j < NNUC; ++j)
						if (val != aux_arr[i][j]) {
							uint64_t bid_new = index_val.bid - val + aux_arr[i][j];
							IndexType index_ret;
							if (!insert(index_dict, bid_new, ValueType(index_val.vid, mask), index_ret)) {
								if (verbose) printf("max_mismatch %d is too high. Reset to %d.\n", max_mismatch, cur_mismatch - 1);
								return;
							}
							buffer2->emplace(IndexType(bid_new, index_val.vid, mask));
						}
				}
			}
			buffer1->pop();
		}
		std::swap(buffer1, buffer2);
		++cur_mismatch;
	}
}


void parse_sample_sheet(const std::string& sample_sheet_file, int& n_barcodes, int& barcode_len, HashType& index_dict, std::vector<std::string>& index_names, int max_mismatch = 1, bool verbose = true) {
	// Assume barcode list does not include Ns
	std::string line;

	n_barcodes = 0;
	barcode_len = 0;
	index_dict.clear();
	index_names.clear();

	bool is_first_line = true;
	std::queue<IndexType>* buffer1 = new std::queue<IndexType>();
	std::queue<IndexType>* buffer2 = new std::queue<IndexType>();

	if (sample_sheet_file.length() > 3 && sample_sheet_file.substr(sample_sheet_file.length() - 3, 3) == ".gz") { // input sample sheet is gzipped
		iGZipFile gin(sample_sheet_file);
		while (gin.next(line)) {
			if (is_first_line) {
				skip_bom(line);
				is_first_line = false;
			}
			parse_one_line(line, n_barcodes, barcode_len, index_dict, index_names, max_mismatch, buffer1);
		}
	}
	else {
		std::ifstream fin(sample_sheet_file);
		while (std::getline(fin, line)) {
			if (is_first_line) {
				skip_bom(line);
				is_first_line = false;
			}
			parse_one_line(line, n_barcodes, barcode_len, index_dict, index_names, max_mismatch, buffer1);
		}
		fin.close();
	}

	insert_mutated_index(index_dict, barcode_len, max_mismatch, buffer1, buffer2, verbose);

	if (verbose) printf("%s is parsed. n_barcodes = %d, and barcode_len = %d.\n", sample_sheet_file.c_str(), n_barcodes, barcode_len);

	int n_amb = 0;
	for (auto&& kv : index_dict)
		if (kv.second.vid < 0) ++n_amb;
	if (verbose) printf("In the index, %d out of %d items are ambigious, percentage = %.2f%%.\n", n_amb, (int)index_dict.size(), n_amb * 100.0 / index_dict.size());
}

#endif
