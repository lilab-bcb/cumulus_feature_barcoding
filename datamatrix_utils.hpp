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
#include <H5Cpp.h>

#include "barcode_utils.hpp"
#include "umi_utils.hpp"

typedef std::unordered_map<int, UMI2Count> Feature2UMI;
typedef std::unordered_map<int, Feature2UMI> Cell2Feature;

const hsize_t CHUNK_SIZE = 80000;

// Either data or (data_size + fillvalue) are specified
void _create_h5_string_dataset(H5::Group& group, const std::string& name, const std::vector<std::string>* data, size_t str_len, size_t data_size, std::string fillvalue) {
	H5::DSetCreatPropList prop_list;
	hsize_t dims_attr[1];

	if (data == nullptr) dims_attr[0] = data_size;
	else dims_attr[0] = data->size();

	H5::DataSpace dataspace_attr(1, dims_attr);
	H5::StrType attr_type(H5::PredType::C_S1, str_len);

	hsize_t chunk_dims[1] = {std::min(CHUNK_SIZE, dims_attr[0])};
	prop_list.setChunk(1, chunk_dims);
	prop_list.setDeflate(6);  // Compression level (0-9)

	if (data == nullptr)
		prop_list.setFillValue(attr_type, fillvalue.c_str());

	H5::DataSet dataset_attr = group.createDataSet(name, attr_type, dataspace_attr, prop_list);

	if (data != nullptr) {
		char* buffer = new char[data->size() * str_len];
		memset(buffer, 0, data->size() * str_len);

		for (size_t i = 0; i < data->size(); ++i)
			strncpy(buffer + i * str_len, (*data)[i].c_str(), (*data)[i].length());

		dataset_attr.write(buffer, attr_type);
		delete[] buffer;
	}
}

void _create_h5_uint_dataset(H5::Group& group, const std::string& name, const std::vector<int>& data) {
	H5::DSetCreatPropList prop_list;
	hsize_t dims_data[] = {data.size()};
	H5::DataSpace dataspace_data(1, dims_data);

	hsize_t chunk_dims[1] = {std::min(CHUNK_SIZE, dims_data[0])};
	prop_list.setChunk(1, chunk_dims);
	prop_list.setDeflate(6);  // Compression level (0-9)
	prop_list.setShuffle();  // Enable shuffle filter (better compression for numeric data)

	H5::DataSet dataset_data = group.createDataSet(name, H5::PredType::NATIVE_UINT, dataspace_data, prop_list);
	dataset_data.write(data.data(), H5::PredType::NATIVE_UINT);
}

void write_molecule_info_h5(
	const std::string& output_name,
	const std::vector<int>& barcode_idx,
	const std::vector<std::string>& barcodes,
	const std::vector<int>& feature_idx,
	const std::vector<std::string>& features,
	int feature_str_len,
	const std::vector<std::string>& umi_names,
	std::vector<int>& umi_counts,
	int umi_len
) {
	const H5std_string filename(output_name + ".molecule_info.h5");

	try {
		H5::H5File fout(filename, H5F_ACC_TRUNC);

		H5::Group grp = fout.openGroup("/");

		// barcodes
		_create_h5_uint_dataset(grp, "barcode_idx", barcode_idx);
		_create_h5_string_dataset(grp, "barcodes", &barcodes, barcodes[0].length(), 0, "");

		// features
		_create_h5_uint_dataset(grp, "feature_idx", feature_idx);
		_create_h5_string_dataset(grp, "features", &features, feature_str_len, 0, "");

		// umis
		_create_h5_string_dataset(grp, "umi", &umi_names, umi_len, 0, "");
		_create_h5_uint_dataset(grp, "count", umi_counts);

	} catch (H5::Exception& error) {
		error.printErrorStack();
		exit(-1);
	}

	printf("%s is written.\n", filename.c_str());
}

void write_count_matrix_h5(
	const std::string& output_name,
	const std::vector<int>& csr_data,
	const std::vector<int>& csr_indices,
	const std::vector<int>& csr_indptr,
	int total_cells,
	int total_features,
	const std::vector<std::string>& barcodes,
	const std::vector<std::string>& features,
	int feature_str_len,
	const std::string& feature_type
) {
	const H5std_string filename(output_name + ".h5");

	try {
		H5::H5File fout(filename, H5F_ACC_TRUNC);

		H5::Group grp = fout.createGroup("/matrix");

		// shape
		hsize_t dims_shape[] = {2};
		H5::DataSpace dataspace_shape(1, dims_shape);
		H5::DataSet dataset_shape = grp.createDataSet("shape", H5::PredType::STD_I32LE, dataspace_shape);
		int shape[2] = {total_features, total_cells};
		dataset_shape.write(shape, H5::PredType::NATIVE_INT);

		// barcodes
		_create_h5_string_dataset(grp, "barcodes", &barcodes, barcodes[0].length(), 0, "");

		// count matrix
		_create_h5_uint_dataset(grp, "data", csr_data);
		_create_h5_uint_dataset(grp, "indices", csr_indices);
		_create_h5_uint_dataset(grp, "indptr", csr_indptr);

		// Features
		H5::Group feature_grp = grp.createGroup("features");

		// name
		_create_h5_string_dataset(feature_grp, "name", &features, feature_str_len, 0, "");
		// id
		_create_h5_string_dataset(feature_grp, "id", &features, feature_str_len, 0, "");

		// feature_type
		std::string ftype = feature_type == "antibody" ? "Antibody Capture" : "CRISPR Guide Capture";
		_create_h5_string_dataset(feature_grp, "feature_type", nullptr, ftype.length(), features.size(), ftype);

		// genome
		_create_h5_string_dataset(feature_grp, "genome", nullptr, 1, features.size(), "");
		_create_h5_string_dataset(feature_grp, "_all_tag_keys", nullptr, 6, 1, "genome");

	} catch (H5::Exception& error) {
		error.printErrorStack();
		exit(-1);
	}

	printf("%s is written.\n", filename.c_str());
}

class DataCollector {
public:
	DataCollector() { clear(); }

	void clear() { data_container.clear(); }

	void insert(int cell_id, int feature_id, uint64_t umi) {
		++data_container[cell_id][feature_id][umi];
	}

	void output(
		const std::string& output_name,
		const std::string& feature_type,
		int feature_start,
		int feature_end,
		const std::vector<std::string>& cell_names,
		int umi_len,
		const std::vector<std::string>& feature_names,
		std::ofstream& freport,
		int n_threads,
		bool verbose_report = true,
		bool is_raw = true
	) {
		std::vector<int> cell_ids;
		std::ofstream fout;

		int total_cells = 0, total_reads = 0, total_umis = 0;
		int total_features = feature_end - feature_start;

		cell_ids.clear();
		for (auto&& kv : data_container)
			cell_ids.push_back(kv.first);
		total_cells = cell_ids.size();
		if (total_cells > 1) std::sort(cell_ids.begin(), cell_ids.end());

		std::vector<int> dummy(total_cells, 0), tot_umis(total_cells, 0);
		std::vector<std::vector<int> > ADTs(total_features, dummy);

		// Collect molecule info and stats
		std::vector<int> barcode_idx;
		std::vector<int> feature_idx;
		std::vector<std::string> umi_names;
		std::vector<int> umi_counts;
		for (int i = 0; i < total_cells; ++i) {
			auto& one_cell = data_container[cell_ids[i]];
			for (auto&& kv1 : one_cell) {
				for (auto&& kv2 : kv1.second) {
					barcode_idx.push_back(i);
					feature_idx.push_back(kv1.first - feature_start);
					umi_names.push_back(binary_to_barcode(kv2.first, umi_len));
					umi_counts.push_back(kv2.second);
					total_reads += kv2.second;
					++total_umis;
					++ADTs[kv1.first - feature_start][i];
					++tot_umis[i];
				}
			}
		}

		// Create a CSR sparse matrix.
		std::vector<int> csr_data;
		std::vector<int> csr_indices;
		std::vector<int> csr_indptr;
		std::vector<std::string> barcodes;

		int indptr = 0;
		csr_indptr.push_back(indptr);
		for (int i = 0; i < total_cells; ++i) {
			if (tot_umis[i] == 0) continue;
			barcodes.push_back(cell_names[cell_ids[i]]);
			for (int j = feature_start; j < feature_end; ++j)
				if (ADTs[j - feature_start][i] > 0) {
					csr_data.push_back(ADTs[j - feature_start][i]);
					csr_indices.push_back(j - feature_start);
					++indptr;
				}
			csr_indptr.push_back(indptr);
		}

		// Collect feature info
		std::vector<std::string> features;
		size_t max_feature_name_len = 0;
		for (int i = feature_start; i < feature_end; ++i) {
			features.push_back(feature_names[i]);
			if (feature_names[i].length() > max_feature_name_len)
				max_feature_name_len = feature_names[i].length();
		}

		write_molecule_info_h5(output_name, barcode_idx, barcodes, feature_idx, features, max_feature_name_len, umi_names, umi_counts, umi_len);
		write_count_matrix_h5(output_name, csr_data, csr_indices, csr_indptr, total_cells, total_features, barcodes, features, max_feature_name_len, feature_type);

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

	void correct_umi(int umi_length, std::string& method) {
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
