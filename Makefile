.PHONY : all clean

all : generate_count_matrix_ADTs

generate_count_matrix_ADTs : generate_count_matrix_ADTs.cpp gzip_utils.hpp barcode_utils.hpp datamatrix_utils.hpp external/slw287r_trimadap/igzip_lib.h external/slw287r_trimadap/izlib.h external/kseq.h
	g++ --std=c++11 -O3 -Wall $< -o $@ -lisal

clean :
	rm -f generate_count_matrix_ADTs
