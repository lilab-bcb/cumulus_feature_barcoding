.PHONY : all clean

all : generate_count_matrix_ADTs

generate_count_matrix_ADTs : generate_count_matrix_ADTs.cpp gzip_utils.hpp barcode_utils.hpp datamatrix_utils.hpp ReadParser.hpp compress.hpp external/slw287r_trimadap/izlib.h external/kseq.h external/concurrentqueue.h external/thread_utils.hpp
	g++ --std=c++11 -O3 -Wall $< -o $@ -lisal -ldeflate -lpthread

clean :
	rm -f generate_count_matrix_ADTs
