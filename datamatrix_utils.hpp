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

typedef std::unordered_map<uint64_t, int> UMI2Count;
typedef std::unordered_map<int, UMI2Count> Feature2UMI;
typedef std::unordered_map<int, Feature2UMI> Cell2Feature;

typedef std::unordered_map<uint64_t, int> UMITable;
typedef UMITable::iterator UMITableIter;

class DisjointSet {
	private:
		std::vector<uint64_t> bid;
		std::vector<int> count;
		std::vector<int> parent;
		std::vector<int> rank;
		std::vector<int> max_count;

		std::vector<int> path;
	public:
		DisjointSet() {}

		void init(UMI2Count umi2count) {
			int n = umi2count.size();

			bid = std::vector<uint64_t>(n, 0);
			count = std::vector<int>(n, 0);
			parent = std::vector<int>(n, -1);
			rank = std::vector<int>(n, 0);
			max_count = std::vector<int>(n, 0);

			int i = 0;
			for (auto& kv: umi2count) {
				bid[i] = kv.first;
				count[i] = kv.second;
				parent[i] = i;
				max_count[i] = i;
				++i;
			}
		}

		std::vector<uint64_t>& get_bid() {
			return bid;
		}

		int find_set(int i) { // with path compression optimization
			path.clear();
			while (i != parent[i]) {
				path.push_back(i);
				i = parent[i];
			}
			for (int& v: path)
				parent[v] = i;
			return i;
		}

		void union_sets(int a, int b) {
			a = find_set(a);
			b = find_set(b);

			if (a != b) {
				if (rank[a] <= rank[b]) {
					parent[a] = b;
					if ((count[max_count[a]] > count[max_count[b]]) || (count[max_count[a]] == count[max_count[b]] && bid[max_count[a]] < bid[max_count[b]]))
						max_count[b] = max_count[a];
				} else if (rank[a] > rank[b]) {
					parent[b] = a;
					if ((count[max_count[a]] < count[max_count[b]]) || (count[max_count[a]] == count[max_count[b]] && bid[max_count[a]] > bid[max_count[b]]))
						max_count[a] = max_count[b];
				}
				if (rank[a] == rank[b])
					++rank[b];
			}
		}

		bool connected(int a, int b) {
			return find_set(a) == find_set(b);
		}

		UMI2Count get_roots() {
			UMI2Count res;
			UMI2Count::iterator it;

			for (int i = 0; i < bid.size(); ++i) {
				int root = find_set(i);
				it = res.find(bid[max_count[root]]);
				if (it != res.end())
					it->second += count[i];
				else
					res.insert(std::make_pair(bid[max_count[root]], count[i]));
			}
			return res;
		}

		void clear() {
			bid.clear();
			count.clear();
			parent.clear();
			rank.clear();
			max_count.clear();
		}
};

class UMICorrectSet {
	private:
		std::vector<UMITable> match_tables;
		DisjointSet ds;
		int umi_length;

		void insert_or_merge() {
			int n_umis = this->ds.get_bid().size();
			for (int i = 0; i < n_umis; ++i)
				for (int j = 0; j < this->umi_length; ++j) {
					uint64_t bid_sub = this->ds.get_bid()[i] & (~aux_arr[j][NNUC]);
					std::pair<UMITableIter, bool> ret;
					ret = match_tables[j].insert(std::make_pair(bid_sub, i));
					if (!ret.second)
						this->ds.union_sets(ret.first->second, i);
				}
		}

	public:
		UMICorrectSet(int umi_length): umi_length(umi_length) {
			this->match_tables = std::vector<UMITable>(umi_length, UMITable());
		}

		void clear() {
			for (int i = 0; i < this->umi_length; ++i)
				this->match_tables[i].clear();
			this->ds.clear();
		}

		void build_set(UMI2Count umi2count) {
			this->ds.init(umi2count);
			insert_or_merge();
		}

		UMI2Count get_corrected_umi_counts() {
			return this->ds.get_roots();
		}

};

class DataCollector {
public:
	DataCollector() { clear(); }

	void clear() { data_container.clear(); }

	void insert(int cell_id, int feature_id, uint64_t umi) {
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
		gout.write("Barcode,Feature,UMI,Count\n");
		for (int i = 0; i < total_cells; ++i) {
			auto& one_cell = data_container[cell_ids[i]];
			for (auto&& kv1 : one_cell) {
				for (auto&& kv2 : kv1.second) {
					gout.write(cell_names[cell_ids[i]]); gout.write(',');
					gout.write(feature_names[kv1.first]); gout.write(',');
					gout.write(binary_to_barcode(kv2.first, umi_len)); gout.write(',');
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

	void correct_umi(int umi_length, const std::vector<std::string>& cell_names, const std::vector<std::string>& feature_names) {
		UMICorrectSet ucs(umi_length);
		for (auto& p: data_container) {
			for (auto& kv: p.second) {
				if (kv.second.size() == 1)
					continue;
				ucs.clear();
				ucs.build_set(kv.second);
				data_container[p.first][kv.first] = ucs.get_corrected_umi_counts();
			}
		}
	}

private:
	Cell2Feature data_container;
};

#endif
