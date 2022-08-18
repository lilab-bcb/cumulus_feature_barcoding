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
	int status; // 1, can proceed; 0, read sequence length 0; -1 EOF; -2 truncated quality string; -3, error reading stream
	gzFile sequence_file = nullptr;
	kseq_t* sequence_kseq = nullptr;


	iGZipFile(const std::string& input_file) : input_file_(input_file), status(1) {
		sequence_file = gzopen(input_file_.c_str(), "r");
		if (sequence_file == nullptr) {
			printf("Cannot find sequence file %s!\n", input_file_.c_str());
			exit(-1);
		}
		sequence_kseq = kseq_init(sequence_file);
	}

	iGZipFile(iGZipFile&& o) {
		input_file_ = std::move(o.input_file_);
		status = o.status;
		sequence_file = o.sequence_file;
		sequence_kseq = o.sequence_kseq;
		o.sequence_file = nullptr;
		o.sequence_kseq = nullptr;
	}

	~iGZipFile() {
		if (sequence_kseq != nullptr) kseq_destroy(sequence_kseq);
		if (sequence_file != nullptr) gzclose(sequence_file);
	}


	std::string get_input_file() const { return input_file_; }

	bool eof() const { return status == -1; }

	bool next(Read& aread) {
		if (status != 1) return false;

		int length = kseq_read(sequence_kseq);

		if (length > 0) {
			aread.name = sequence_kseq->name.s;
			aread.comment = sequence_kseq->comment.s;
			aread.seq = sequence_kseq->seq.s;
			aread.qual = sequence_kseq->is_fastq ? sequence_kseq->qual.s : "";
			return true;
		}

		status = length;
		return false;
	}


	bool next(std::string& line) { // Only get one line
		if (status != 1) return false;

		int ret = ks_getuntil(sequence_kseq->f, KS_SEP_LINE, &sequence_kseq->comment, 0);

		if (ret >= 0) {
			line = sequence_kseq->comment.s;
			return true;
		}

		status = ret;
		return false;
	}


	bool check_error(int cnt, bool is_read = true) const { // cnt: 0-based count; return true if there is an error.
		switch(status) {
			case -3:
				printf("File %s, %s %d (numbered from 1): error reading stream, file might be corrupted!\n", input_file_.c_str(), (is_read ? "read" : "line"), cnt + 1);
				return true;
			case -2:
				printf("File %s, %s %d (numbered from 1): truncated quality string!\n", input_file_.c_str(), (is_read ? "read" : "line"), cnt + 1);
				return true;
			case 0:
				printf("File %s, %s %d (numbered from 1): Detected a read with 0 sequence length!\n", input_file_.c_str(), (is_read ? "read" : "line"), cnt + 1);
				return true;
		}
		return false;
	}
};


struct oGZipFile {
	FILE *fo = nullptr;
	Compressor *compressor = nullptr;

	oGZipFile(const std::string& output_file, int num_threads = 1, bool bgzf = false, size_t buffer_size = compressor_buffer_size, int compression_level = 6) {
		fo = fopen(output_file.c_str(), "wb");
		if (fo == nullptr) {
			printf("Cannot creat output file %s!\n", output_file.c_str());
			exit(-1);
		}

		assert(num_threads >= 1);
		compressor = nullptr;
		if (num_threads == 1) {
			compressor = new SingleThreadCompressor(buffer_size, compression_level, bgzf);
		} 
		else {
			compressor = new MultiThreadsCompressor(num_threads, buffer_size, compression_level, bgzf);
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
			flush(true);
			fclose(fo);
			delete compressor;
			fo = nullptr;
			compressor = nullptr;			
		}
	}

	void flush(bool last_flush = false) {
		size_t cprs_size = compressor->compress();
		compressor->flushOut(fo, cprs_size, last_flush);
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
