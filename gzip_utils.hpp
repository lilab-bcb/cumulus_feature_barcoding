#ifndef GZIP_UTILS
#define GZIP_UTILS

#include <cctype>
#include <cassert>
#include <string>
#include <cstdio>

#include "external/slw287r_trimadap/izlib.h"
#include "external/kseq.h"
#include "compress.hpp"


struct Read {
	std::string name, comment, seq, qual; // qual.empty() -> FASTA reads

	size_t size() const {
		return name.length() + (comment.empty() ? 0 : 1 + comment.length()) + seq.length() + 3 + (qual.empty() ? 0 : 3 + qual.length());
	}

	std::string toString() const {
		return "@" + name + (comment.empty() ? "" : (" " + comment)) + "\n" + seq + "\n" + (qual.empty() ? "" : "+\n" + qual + "\n");
	}
};

KSEQ_INIT(gzFile, gzread);

struct iGZipFile {
	std::string input_file_;
	gzFile sequence_file = nullptr;
	kseq_t* sequence_kseq = nullptr;


	iGZipFile(const std::string& input_file) : input_file_(input_file) {
		sequence_file = gzopen(input_file_.c_str(), "r");
		if (sequence_file == nullptr) {
			printf("Cannot find sequence file %s!\n", input_file_.c_str());
			exit(-1);
		}
		sequence_kseq = kseq_init(sequence_file);
	}

	iGZipFile(iGZipFile&& o) {
		input_file_ = std::move(o.input_file_);
		sequence_file = o.sequence_file;
		sequence_kseq = o.sequence_kseq;
		o.sequence_file = nullptr;
		o.sequence_kseq = nullptr;
	}

	~iGZipFile() {
		if (sequence_kseq != nullptr) kseq_destroy(sequence_kseq);
		if (sequence_file != nullptr) gzclose(sequence_file);
	}


	bool next(Read& aread) {
		int length = kseq_read(sequence_kseq);
		if (length > 0) {
			aread.name = sequence_kseq->name.s;
			aread.comment = sequence_kseq->comment.s;
			aread.seq = sequence_kseq->seq.s;
			aread.qual = sequence_kseq->is_fastq ? sequence_kseq->qual.s : "";
		}
		else if (length == -1) 
			return false; // End of file
		else if (length == -3) {
			printf("File %s: error reading stream; didn't reach the end of sequence file, which might be corrupted!\n", input_file_.c_str());
			exit(-1);			
		}
		else if (length == -2) {
			printf("File %s: truncated quality string!\n", input_file_.c_str());
			exit(-1);			
		}
		else {
			assert(length == 0);
			printf("File %s: Detected a read with 0 sequence length!\n", input_file_.c_str());
			exit(-1);
		}

		return true;
	}


	bool next(std::string& line) { // Only get one line
		int ret = ks_getuntil(sequence_kseq->f, KS_SEP_LINE, &sequence_kseq->comment, 0);

		if (ret == -3) {
			printf("File %s: error reading stream in next(line)!\n", input_file_.c_str());
			exit(-1);						
		}

		if (ret >= 0) {
			line = sequence_kseq->comment.s;
			return true;
		}

		return false;
	}
};


struct oGZipFile {
	FILE *fo = nullptr;
	Compressor *compressor = nullptr;

	oGZipFile(const std::string& output_file, int num_threads = 1, size_t buffer_size = compressor_buffer_size, int compression_level = 6) {
		fo = fopen(output_file.c_str(), "wb");
		if (fo == nullptr) {
			printf("Cannot creat output file %s!\n", output_file.c_str());
			exit(-1);
		}

		assert(num_threads >= 1);
		compressor = nullptr;
		if (num_threads == 1) {
			compressor = new SingleThreadCompressor(buffer_size, compression_level);
		} 
		else {
			compressor = new MultiThreadsCompressor(num_threads, buffer_size, compression_level);
		}
	}

	oGZipFile(oGZipFile&& o) {
		fo = o.fo;
		compressor = o.compressor;
		o.fo = nullptr;
		o.compressor = nullptr;
	}

	~oGZipFile() {
		close();
	}

	void close() {
		if (fo != nullptr) {
			flush();
			fclose(fo);
			delete compressor;
			fo = nullptr;
			compressor = nullptr;			
		}
	}

	void flush() {
		size_t cprs_size = compressor->compress();
		if (cprs_size > 0) compressor->flushOut(fo, cprs_size);
	}

	void write(const Read& aread) {
		if (compressor->needFlush(aread.size())) flush();
		compressor->writeToBuffer('@');
		compressor->writeToBuffer(aread.name.c_str(), aread.name.length());
		if (!aread.comment.empty()) {
			compressor->writeToBuffer(' ');
			compressor->writeToBuffer(aread.comment.c_str(), aread.comment.length());
		}
		compressor->writeToBuffer('\n');
		compressor->writeToBuffer(aread.seq.c_str(), aread.seq.length());
		compressor->writeToBuffer('\n');
		if (!aread.qual.empty()) {
			compressor->writeToBuffer('+');
			compressor->writeToBuffer('\n');
			compressor->writeToBuffer(aread.qual.c_str(), aread.qual.length());
			compressor->writeToBuffer('\n');
		}
	}

	void write(const std::string& line) {
		if (compressor->needFlush(line.length())) flush();
		compressor->writeToBuffer(line.c_str(), line.length());
	}

	void write(char c) {
		if (compressor->needFlush(1)) flush();
		compressor->writeToBuffer(c);
	}
};

#endif
