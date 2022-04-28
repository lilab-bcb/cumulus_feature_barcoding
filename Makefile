.PHONY : all clean

all : generate_count_matrix_ADTs

generate_count_matrix_ADTs : generate_count_matrix_ADTs.cpp gzip_utils.hpp barcode_utils.hpp datamatrix_utils.hpp kseq.h
	g++ --std=c++17 -O3 -lz -Wall $< -o $@ FQFeeder/build/src/libfqfeeder.a

clean :
	rm -f generate_count_matrix_ADTs
