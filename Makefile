.PHONY : all clean

all : generate_count_matrix_ADTs

generate_count_matrix_ADTs : generate_count_matrix_ADTs.cpp gzip_utils.hpp barcode_utils.hpp datamatrix_utils.hpp kseq.h
	g++ --std=c++11 -O3 -lz $< -o $@

clean :
	rm -f generate_count_matrix_ADTs
