#ifndef DATAMATRIX_UTILS
#define DATAMATRIX_UTILS

#include <cstdio>
#include <cstdint>
#include <string>
#include <vector>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <unordered_map>

#include "gzip_utils.hpp"
#include "barcode_utils.hpp"
#include "umi_utils.hpp"

typedef std::unordered_map<int, UMI2Count> Feature2UMI;
typedef std::unordered_map<int, Feature2UMI> Cell2Feature;


class DataCollector {
public:
	DataCollector() { clear(); }

	void clear() { data_container.clear(); }

	void insert(int cell_id, int feature_id, uint64_t umi) {
		++data_container[cell_id][feature_id][umi];
	}

	void output(const std::string& output_name, const std::string& feature_type, int feature_start, int feature_end, const std::vector<std::string>& cell_names, int umi_len, const std::vector<std::string>& feature_names, std::ofstream& freport, int n_threads, bool verbose_report = true, bool is_raw = true) {
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

		if (is_raw) {
			report_buffer << std::endl<< "Section "<< output_name<< std::endl;
			report_buffer << "Number of valid cell barcodes: "<< total_cells<< std::endl;
			report_buffer << "Number of valid reads (with matching cell and feature barcodes): "<< total_reads<< std::endl;
			report_buffer << "Mean number of valid reads per cell barcode: "<< std::fixed<< std::setprecision(2)<< (total_cells > 0 ? total_reads * 1.0 / total_cells : 0.0)<< std::endl;
			report_buffer << "Number of valid UMIs (with matching cell and feature barcodes): "<< total_umis<< std::endl;
			report_buffer << "Mean number of valid UMIs per cell barcode: "<< std::fixed<< std::setprecision(2)<< (total_cells > 0 ? total_umis * 1.0 / total_cells : 0.0)<< std::endl;
			report_buffer << "Sequencing saturation: "<< std::fixed<< std::setprecision(2)<< (total_reads > 0 ? 100.0 - total_umis * 100.0 / total_reads : 0.0)<< "%"<< std::endl;
		} else {
			report_buffer << "After UMI correction:" << std::endl;
			report_buffer << "\tNumber of valid UMIs (with matching cell and feature barcodes): "<< total_umis<< std::endl;
			report_buffer << "\tMean number of valid UMIs per cell barcode: "<< std::fixed<< std::setprecision(2)<< (total_cells > 0 ? total_umis * 1.0 / total_cells : 0.0)<< std::endl;
			report_buffer << "\tSequencing saturation: "<< std::fixed<< std::setprecision(2)<< (total_reads > 0 ? 100.0 - total_umis * 100.0 / total_reads : 0.0)<< "%"<< std::endl;
		}

		if (verbose_report) {
			freport << report_buffer.str();
		}
	}

	void correct_umi(int umi_length, const std::vector<std::string>& cell_names, const std::vector<std::string>& feature_names, std::string& method) {
		UMICorrectSet ucs(umi_length);
		for (auto& p: data_container) {
			for (auto& kv: p.second) {
				if (kv.second.size() == 1)
					continue;
				ucs.clear();
				ucs.build_set(kv.second, method);
				data_container[p.first][kv.first] = ucs.get_corrected_umi_counts();
			}
		}
	}

private:
	Cell2Feature data_container;
	std::ostringstream report_buffer;
};

#endif
