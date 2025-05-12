#ifndef AUTO_UTILS
#define AUTO_UTILS

#include <string>
#include <vector>
#include <unordered_map>
#include <filesystem>

#include "barcode_utils.hpp"
#include "ReadParser.hpp"


const int totalseq_A_pos = 0;
const int totalseq_BC_pos = 10;


std::unordered_map<std::string, std::vector<std::string>> compound_chemistry_dict = {
	{"auto", {"10x_v2", "SC3Pv3:Poly-A", "SC3Pv3:CS1", "SC3Pv4:Poly-A", "SC3Pv4:CS1", "SC5Pv3", "multiome"}},
	{"threeprime", {"10x_v2", "SC3Pv3:Poly-A", "SC3Pv3:CS1", "SC3Pv4:Poly-A", "SC3Pv4:CS1"}},
	{"fiveprime", {"10x_v2", "SC5Pv3"}},
	{"SC3Pv3", {"SC3Pv3:Poly-A", "SC3Pv3:CS1"}},
	{"SC3Pv4", {"SC3Pv4:Poly-A", "SC3Pv4:CS1"}},
};

std::unordered_map<std::string, std::string> cb_inclusion_file_dict = {
	{"10x_v2", "737K-august-2016.txt.gz"},
	{"SC3Pv2", "737K-august-2016.txt.gz"},
	{"SC5Pv2", "737K-august-2016.txt.gz"},
	{"SC3Pv3:Poly-A", "3M-february-2018_TRU.txt.gz"},
	{"SC3Pv3:CS1", "3M-february-2018_NXT.txt.gz"},
	{"SC3Pv4:Poly-A", "3M-3pgex-may-2023_TRU.txt.gz"},
	{"SC3Pv4:CS1", "3M-3pgex-may-2023_NXT.txt.gz"},
	{"SC5Pv3", "3M-5pgex-jan-2023.txt.gz"},
	{"multiome", "737K-arc-v1.txt.gz"},
};

inline std::string safe_substr(const std::string& sequence, int pos, int length) {
	if (pos + length > sequence.length()) {
		printf("Error: Sequence length %d is too short (expected to be at least %d)!\n", (int)sequence.length(), pos + length);
		exit(-1);
	}
	return sequence.substr(pos, length);
}

void auto_detection(
    std::string& cb_dir,
    std::vector<std::vector<std::string>>& inputs,
    std::string& feature_type,
    int& feature_blen,
    HashType& feature_index,
    int& umi_len,
    std::string& chemistry,
    std::string& totalseq_type,
    int& max_mismatch_cell,
    int& barcode_pos,
    std::string& scaffold_sequence,
    bool verbose = true
) {
	// Redirect 10x barcode inclusion list files
	for (auto& p : cb_inclusion_file_dict) {
		p.second = cb_dir + p.second;
	}

	const int nskim = 10000;  // Look at first 10,000 reads if auto-detection is needed.
	int cnt;  // Hide the global "cnt" which is for total numbers of reads.
	std::size_t pos;

	auto it = compound_chemistry_dict.find(chemistry);
	if (it == compound_chemistry_dict.end()) {    // The given chemistry name is not a compound chemistry type
		if (cb_inclusion_file_dict.find(chemistry) == cb_inclusion_file_dict.end()) {
			printf("Error: Unknown chemistry type: %s!\n", chemistry.c_str());
			exit(-1);
		} else if (chemistry == "10x_v2") {
			printf("Error: %s chemistry type is for internal use only!\n", chemistry.c_str());
			exit(-1);
		}
	} else {
		// A compound chemistry is given. Auto-detect
		int n_chems = it->second.size();
		std::string cur_chem, cb_filename;
		int n_cb, len_cb;
		std::vector<std::string> dummy;  // Placeholder. Not used.
		uint64_t binary_cb;
		Read read1;

		//Build count map
		std::vector<std::string> chem_names = std::vector<std::string>(n_chems, "");
		std::vector<int> chem_cnts = std::vector<int>(n_chems, 0);
		std::vector<HashType> chem_cb_indexes = std::vector<HashType>(n_chems, HashType());
		for (int i = 0; i < n_chems; ++i) {
			cur_chem = it->second[i];
			chem_names[i] = cur_chem;
            cb_filename = cb_inclusion_file_dict[cur_chem];
            if (verbose)  printf("[Auto-detection] Loading %s cell barcode file...\n", cur_chem.c_str());
            if (!std::filesystem::exists(std::filesystem::path(cb_filename)))
                printf("Warning: Cannot find cell barcode file '%s' for chemistry %s!\n", cb_filename.c_str(), cur_chem.c_str());
			parse_sample_sheet(cb_filename, n_cb, len_cb, chem_cb_indexes[i], dummy, 0, false);
		}

		// Count cell barcode matches
		cnt = 0;
		for (auto&& input_pair : inputs) {
			iGZipFile gzip_in_r1(input_pair[0]);
			while (gzip_in_r1.next(read1) && cnt < nskim) {
				binary_cb = barcode_to_binary(safe_substr(read1.seq, 0, len_cb));
				for (int i = 0; i < n_chems; ++i) {
					if (chem_cb_indexes[i].find(binary_cb) != chem_cb_indexes[i].end())
						++chem_cnts[i];
				}
				++cnt;
			}
			if (cnt == nskim) break;
		}

		// Decide chemistry based on count results
		int max_cnt = -1;
		int snd_max_cnt = -1;
		std::string chem_max;
		std::string chem_snd_max;

		for (int i = 0; i < n_chems; ++i)
			if (chem_cnts[i] > max_cnt) {
				snd_max_cnt = max_cnt;
				chem_snd_max = chem_max;
				max_cnt = chem_cnts[i];
				chem_max = chem_names[i];
			} else if (chem_cnts[i] > snd_max_cnt) {
				snd_max_cnt = chem_cnts[i];
				chem_snd_max = chem_names[i];
			}

		if (max_cnt > 0) {
			if (snd_max_cnt > 0) {
                if (verbose)
				    printf("[Auto-detection] Top 2 chemistry types with maximum matches on cell barcodes in the first %d reads: %s (%d matches), %s (%d matches).\n", nskim, chem_max.c_str(), max_cnt, chem_snd_max.c_str(), snd_max_cnt);
				if (static_cast<float>(max_cnt) / nskim < 0.05) {
					printf("Error: No chemistry has matched reads exceeding 5%% of the first %d reads! Please check if you specify the correct chemistry type, or if it is a 10x assay!\n", nskim);
					exit(-1);
				} else if (static_cast<float>(max_cnt - snd_max_cnt) / nskim < 0.1) {
					printf("Error: Cannot decide chemistry type, as the number of matches of the top 2 chemistry types are too close (< 10%% of the first %d reads)! Please check if it is a 10x assay!\n", nskim);
					exit(-1);
				}
			} else {
                if (verbose)
				    printf("[Auto-detection] Only 1 chemistry type has matches on cell barcodes in the first %d reads: %s (%d matches).\n", nskim, chem_max.c_str(), max_cnt);
            }

			chemistry = chem_max;
		} else {
			printf("Error: Failed at chemistry detection: No cell barcode match in the first %d reads! Please check if it is a 10x assay!", nskim);
			exit(-1);
		}
	} // End of chemistry detection

	// Detect umi_len and max_mismatch_cell
	umi_len = (chemistry == "10x_v2" || chemistry == "SC3Pv2" || chemistry == "SC5Pv2") ? 10 : 12;
	max_mismatch_cell = (chemistry == "10x_v2" || chemistry == "SC3Pv2" || chemistry == "SC5Pv2" || chemistry == "multiome") ? 1 : 0;

	// Detect totalseq_type (for antibody assays) and barcode_pos
	if (feature_type == "hashing" || feature_type == "citeseq" || feature_type == "adt") {
		// Detect totalseq_type
		pos = chemistry.find_first_of(':');
		if (pos != std::string::npos) {
			std::string capture_method = chemistry.substr(pos + 1);
			totalseq_type = capture_method == "Poly-A" ? "TotalSeq-A" : "TotalSeq-B";
		} else if (chemistry == "SC5Pv3")
			totalseq_type = "TotalSeq-C";
		else {
			// 10x_v2 or multiome
			// if specify --barcode-pos, must be a customized assay
			if (barcode_pos < 0) {
				int ntotA, ntotC;
				uint64_t binary_feature;
				Read read2;
				HashIterType feature_iter;

				cnt = ntotA = ntotC = 0;
				for (auto&& input_pair : inputs) {
					iGZipFile gzip_in_r2(input_pair[1]);
					while (gzip_in_r2.next(read2) && cnt < nskim) {
						binary_feature = barcode_to_binary(safe_substr(read2.seq, totalseq_A_pos, feature_blen));
						feature_iter = feature_index.find(binary_feature);
						ntotA += (feature_iter != feature_index.end() && feature_iter->second.vid >= 0);

						if (read2.seq.length() >= totalseq_BC_pos + feature_blen) {
							binary_feature = barcode_to_binary(safe_substr(read2.seq, totalseq_BC_pos, feature_blen));
							feature_iter = feature_index.find(binary_feature);
							ntotC += (feature_iter != feature_index.end() && feature_iter->second.vid >= 0);
						}
						++cnt;
					}
					if (cnt == nskim) break;
				}

				if (verbose)  printf("[Auto-detection] ntotA = %d, ntotC = %d.\n", ntotA, ntotC);
				if (ntotA < 10 && ntotC < 10) {
					printf("Error: Detected less than 10 feature barcodes in the first %d reads! Maybe you should consider to reverse complement your barcodes?\n", nskim);
					exit(-1);
				}
				totalseq_type = (ntotA > ntotC ? "TotalSeq-A" : "TotalSeq-C");
				if (chemistry == "10x_v2")
					chemistry = totalseq_type == "TotalSeq-A" ? "SC3Pv2" : "SC5Pv2";
			}
		}

		// Decide barcode_pos based on TotalSeq type
		barcode_pos = totalseq_type == "TotalSeq-A" ? totalseq_A_pos : totalseq_BC_pos;

	} else if (feature_type == "cmo")  barcode_pos = 0;
    else {
		if (feature_type != "crispr") {
			printf("Error: Do not support unknown feature type %s!\n", feature_type.c_str());
			exit(-1);
		}

		if (barcode_pos < 0 && scaffold_sequence == "") {
			barcode_pos = 0;  // default is 0
            if (verbose)
			    printf("[Auto-detection] Warning: Automatically set barcode start position to %d, as neither --barcode-pos nor --scaffold-sequence is specified.\n", barcode_pos);
		}
	}

    printf("[Auto-detection] Chemistry type detected: %s.\n", chemistry.c_str());
    printf("[Auto-detection] Set UMI length to %d, and set maximum cell barcode mismatch to %d.\n", umi_len, max_mismatch_cell);
    if (feature_type == "hashing" || feature_type == "citeseq" || feature_type == "adt")
        printf("[Auto-detection] TotalSeq type detected: %s. Feature barcodes start from 0-based position %d.\n", totalseq_type.c_str(), barcode_pos);
}

#endif
