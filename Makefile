.PHONY : all clean

UNAME_S = $(shell uname -s)

ifeq ($(UNAME_S), Linux)
	OPT_INCLUDE = /usr/include/hdf5/serial
	OPT_LIB = /usr/lib/x86_64-linux-gnu/hdf5/serial
	CXX = g++
else ifeq ($(UNAME_S), Darwin)
	OPT_INCLUDE = /opt/homebrew/include
	OPT_LIB = /opt/homebrew/lib
	CXX = clang++
else
	$(error Unsupported OS: $(UNAME_S).)
endif

all : generate_feature_barcode_matrices

generate_feature_barcode_matrices : generate_feature_barcode_matrices.cpp gzip_utils.hpp barcode_utils.hpp datamatrix_utils.hpp umi_utils.hpp auto_utils.hpp ReadParser.hpp compress.hpp external/slw287r_trimadap/izlib.h external/kseq.h external/concurrentqueue.h external/thread_utils.hpp
	$(CXX) --std=c++20 -O3 -Wall $< -o $@ -lisal -ldeflate -lpthread -lhdf5 -lhdf5_cpp -L$(OPT_LIB) -I$(OPT_INCLUDE)

clean :
	rm -f generate_feature_barcode_matrices
