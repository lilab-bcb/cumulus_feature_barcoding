#ifndef DATAMATRIX_UTILS
#define DATAMATRIX_UTILS

#include <cstdio>
#include <cstdint>
#include <string>
#include <vector>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <unordered_map>

#include "gzip_utils.hpp"
#include "barcode_utils.hpp"

typedef std::unordered_map<BinaryCodeType, int, BinaryCodeTypeHash> UMI2Count;
typedef std::unordered_map<int, UMI2Count> Feature2UMI;
typedef std::unordered_map<int, Feature2UMI> Cell2Feature;



/*
const int MAX_UMI_LEN = 12;

class DisjointSet {
	private:
		std::unordered_map<uint64_t, uint64_t> parent;
		std::unordered_map<uint64_t, int> count;
		int mismatch_pos;

	public:
		// Constructor: Initializes each element as its own parent (a separate set).
		DisjointSet(int pos) {
			assert (pos >= 0 && pos < MAX_UMI_LEN);
			parent.clear();
			count.clear();
			mismatch_pos = pos;
		}

		int find(int i) {
			if (parent[i] != i) {
				// Path compression: Make the parent of i the root of the set.
				parent[i] = find(parent[i]);
			}
			return parent[i];
		}

		void unite(int x, int y) {
			int rootX = find(x);
			int rootY = find(y);

			if (rootX != rootY) {
				// Attach the shorter tree to the taller tree.
				if (rank[rootX] <= rank[rootY]) {
					parent[rootX] = rootY;
				} else if (rank[rootX] > rank[rootY]) {
					parent[rootY] = rootX;
				}
			}
		}

		bool connected(int x, int y){
			return find(x) == find(y);
		}

		void build_set(std::vector<std::pair<uint64_t, int> > umi_counts) {
			for (auto& p: umi_counts) {
				parent[p.first] = p.first;
				count[p.first] = p.second;
			}

			for (auto& kv: parent) {
				uint64_t val = kv.first & aux_arr[mismatch_pos][NNUC];
				for (int j = 0; j < NNUC; ++j)
					if (val != aux_arr[mismatch_pos][j]) {
						uint64_t bid_new = kv.first - val + aux_arr[mismatch_pos][j];
						if (parent.find(bid_new) != parent.end())
							unite(kv.first, bid_new);
					}
			}
		}

		std::vector<std::pair<uint64_t, int> > get_corrected_umi_counts() {

		}
};
*/


class DataCollector {
public:
	DataCollector() { clear(); }

	void clear() { data_container.clear(); }

	void insert(int cell_id, int feature_id, BinaryCodeType umi) {
		++data_container[cell_id][feature_id][umi];
	}

	void output(const std::string& output_name, const std::string& feature_type, int feature_start, int feature_end, const std::vector<std::string>& cell_names, int umi_len, const std::vector<std::string>& feature_names, std::ofstream& freport, int n_threads) {
		std::vector<int> cell_ids;
		std::ofstream fout;

		int total_cells = 0, total_reads = 0, total_umis = 0;

		cell_ids.clear();
		for (auto&& kv : data_container)
			cell_ids.push_back(kv.first);
		total_cells = cell_ids.size();
		if (total_cells > 1) std::sort(cell_ids.begin(), cell_ids.end());

		std::vector<int> dummy(total_cells, 0), tot_umis(total_cells, 0);
		std::vector<std::vector<int> > ADTs(feature_end - feature_start, dummy);

		oGZipFile gout(output_name + ".stat.csv.gz", n_threads);
		gout.write("Barcode,UMI,Feature,Count\n");
		for (int i = 0; i < total_cells; ++i) {
			auto& one_cell = data_container[cell_ids[i]];
			for (auto&& kv1 : one_cell) {
				for (auto&& kv2 : kv1.second) {
					gout.write(cell_names[cell_ids[i]]); gout.write(',');
					gout.write(binary_to_barcode(kv2.first, umi_len)); gout.write(',');
					gout.write(feature_names[kv1.first]); gout.write(',');
					gout.write(std::to_string(kv2.second)); gout.write('\n');
					total_reads += kv2.second;
					++total_umis;
					++ADTs[kv1.first - feature_start][i];
					++tot_umis[i];
				}
			}
		}
		gout.close();
		printf("%s.stat.csv.gz is written.\n", output_name.c_str());

		fout.open(output_name + ".csv");
		fout<< (feature_type == "antibody" ? "Antibody" : "CRISPR");
		for (int i = 0; i < total_cells; ++i)
			if (tot_umis[i] > 0) fout<< ","<< cell_names[cell_ids[i]];
		fout<< std::endl;
		for (int i = feature_start; i < feature_end; ++i) {
			fout<< feature_names[i];
			for (int j = 0; j < total_cells; ++j)
				if (tot_umis[j] > 0) fout<< ","<< ADTs[i - feature_start][j];
			fout<< std::endl;
		}
		fout.close();
		printf("%s.csv is written.\n", output_name.c_str());

		freport<< std::endl<< "Section "<< output_name<< std::endl;
		freport<< "Number of valid cell barcodes: "<< total_cells<< std::endl;
		freport<< "Number of valid reads (with matching cell and feature barcodes): "<< total_reads<< std::endl;
		freport<< "Mean number of valid reads per cell barcode: "<< std::fixed<< std::setprecision(2)<< (total_cells > 0 ? total_reads * 1.0 / total_cells : 0.0)<< std::endl;
		freport<< "Number of valid UMIs (with matching cell and feature barcodes): "<< total_umis<< std::endl;
		freport<< "Mean number of valid UMIs per cell barcode: "<< std::fixed<< std::setprecision(2)<< (total_cells > 0 ? total_umis * 1.0 / total_cells : 0.0)<< std::endl;
		freport<< "Sequencing saturation: "<< std::fixed<< std::setprecision(2)<< (total_reads > 0 ? 100.0 - total_umis * 100.0 / total_reads : 0.0)<< "%"<< std::endl;
	}
/*
	void correct_umi_counts(int umi_len, int umi_mismatch) {
		std::vector<DisjointSet> ds (umi_len, DisjointSet(int(pow(2, umi_len - 1))));
		std::unordered_map<int, std::vector<std::pair<uint64_t, int> > > feature2umi;

		uint64_t umi_conv;
		for (auto& p: data_container) {
			uint64_t& cell_bid = p.first;
			feature2umi.clear();
			for (auto& p1: p.second) {
				uint64_t& umi_bid = p1.first;
				for (auto& p2: p1.second) {
					int& feature_id = p2.first;
					if (feature2umi.find(feature_id) != feature2umi.end())
						feature2umi[feature_id].emplace_back(std::make_pair(umi_bid, p2.second));
					else
						feature2umi[feature_id] = vector(1, std::make_pair(umi_bid, p2.second));
				}
			}
			for (auto& kv: feature2umi) {
				for (auto& v: kv) {
					for (int i = 0; i < MAX_UMI_LEN; ++i) {
						uint64_t val = v.first & aux_arr[i][NNUC];
						for (int j = 0; j < NNUC; ++j) {
							if (val != aux_arr[i][j]) {
								uint64_t bid_new = v.first - val + aux_arr[i][j];
								ds[i].set_rank()
							}
						}
					}
				}
			}
		}
	}
*/
private:
	Cell2Feature data_container;
};

#endif
