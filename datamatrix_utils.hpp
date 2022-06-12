#ifndef DATAMATRIX_UTILS
#define DATAMATRIX_UTILS

#include <ctime>
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

typedef std::unordered_map<int, int> Feature2Count;
typedef std::unordered_map<uint64_t, Feature2Count> UMI2Feature;
typedef std::unordered_map<int, UMI2Feature> Cell2UMI;

class DataCollector {
public:
	DataCollector() { clear(); }

	void clear() { data_container.clear(); }

	void insert(int cell_id, uint64_t umi, int feature_id) {
		++data_container[cell_id][umi][feature_id];
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

		time_t start_, end_;
		start_ = time(NULL);

		oGZipFile gout(output_name + ".stat.csv.gz", n_threads);
		gout.write("Barcode,UMI,Feature,Count\n");
		for (int i = 0; i < total_cells; ++i) {
			auto& one_cell = data_container[cell_ids[i]];
			for (auto&& kv1 : one_cell) {
				for (auto&& kv2 : kv1.second) {
					gout.write(cell_names[cell_ids[i]]); gout.write(',');
					gout.write(binary_to_barcode(kv1.first, umi_len)); gout.write(',');
					gout.write(feature_names[kv2.first]); gout.write(',');
					gout.write(std::to_string(kv2.second)); gout.write('\n');
					total_reads += kv2.second;
					++total_umis;
					++ADTs[kv2.first - feature_start][i];
					++tot_umis[i];
				}
			}
		}
		gout.close();
		printf("%s.stat.csv.gz is written.\n", output_name.c_str());

		end_ = time(NULL);
		printf("Time spent compressing = %.2fs.\n", difftime(end_, start_));

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

private:
	Cell2UMI data_container;
};

#endif
