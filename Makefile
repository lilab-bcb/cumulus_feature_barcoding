.PHONY : all clean

all : generate_feature_barcode_matrices

generate_feature_barcode_matrices : generate_feature_barcode_matrices.cpp gzip_utils.hpp barcode_utils.hpp datamatrix_utils.hpp umi_utils.hpp auto_utils.hpp ReadParser.hpp compress.hpp external/slw287r_trimadap/izlib.h external/kseq.h external/concurrentqueue.h external/thread_utils.hpp
	g++ --std=c++20 -O3 -Wall $< -o $@ -lisal -ldeflate -lpthread -lhdf5 -lhdf5_cpp

clean :
	rm -f generate_feature_barcode_matrices
