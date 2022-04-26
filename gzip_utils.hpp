#ifndef GZIP_UTILS
#define GZIP_UTILS

#include <cctype>
#include <string>
#include <istream>
#include <streambuf>
#include <fstream>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include <zlib.h>
#include "kseq.h"

struct Read {
	std::string name, seq, qual;

	std::string toString() const {
		return "@" + name + "\n" + seq + "\n+\n" + qual + "\n";
	}
};

struct iGZipKseqFile {
	KSEQ_INIT(gzFile, gzread);

    gzFile sequence_file = NULL;
	kseq_t *sequence_kseq = NULL;
    
	iGZipKseqFile(const std::string& input_file = "") {
		if (input_file != "") open(input_file);
	}

	~iGZipKseqFile() { close(); }

	bool open(const std::string& input_file) {
		sequence_file = gzopen(input_file.c_str(), "r");
		if (sequence_file == NULL) {
			printf("Cannot find sequence file %s!\n", input_file.c_str());
			exit(-1);
		}
		sequence_kseq = kseq_init(sequence_file);
		return true;
	}

	bool close() {
		if (sequence_kseq != NULL) {
			kseq_destroy(sequence_kseq);
			gzclose(sequence_file);
			sequence_kseq = NULL;
		}
		return true;
	}

	int next(Read& aread) {
		int length = kseq_read(sequence_kseq);
		while (length == 0) {  // Skip the sequences of length 0
			length = kseq_read(sequence_kseq);
		}
		if (length > 0) {
			aread.name = sequence_kseq->name.s;
			aread.seq = sequence_kseq->seq.s;
			aread.qual = sequence_kseq->qual.s;
			return 4;
		}
		else {
			if (length != -1) { // make sure to reach the end of file rather than meet an error
				printf("Didn't reach the end of sequence file, which might be corrupted!\n");
				exit(-1);
			}
			return -1;
		}
	}
};

struct iGZipFile {
	std::ifstream fin;
	boost::iostreams::filtering_istream gin;
	std::string line;

	iGZipFile(const std::string& input_file = "") {
		if (input_file != "") open(input_file);
	}

	~iGZipFile() { close(); }

	bool open(const std::string& input_file) {
		size_t pos = input_file.find(".gz");
		bool isGZ = (pos != std::string::npos) && (pos + 3 == input_file.length());

		fin.open(input_file, std::ios_base::in | std::ios_base::binary);
		if (isGZ) gin.push(boost::iostreams::gzip_decompressor());
		gin.push(fin);

		return true;
	}

	bool close() {
		if (!gin.empty()) gin.reset();
		if (fin.is_open()) fin.close();
		return true;
	}

	int next(Read& aread) {
		bool success;

		success = (bool)std::getline(gin, aread.name);
		if (!success) return 0;
		aread.name = aread.name.substr(1);
		success = (bool)std::getline(gin, aread.seq);
		if (!success) return 1;
		success = (bool)std::getline(gin, line);
		if (!success) return 2;
		success = (bool)std::getline(gin, aread.qual);
		return (success ? 4 : 3);
	}

	boost::iostreams::filtering_istream& operator()() {
		return gin;
	}

	bool getline(std::string& line) {
		line.clear();

		std::istream::sentry se(gin, true);
		std::streambuf* sb = gin.rdbuf();

		int c;
		do {
			c = sb->sbumpc();
			if (c == '\n') return true;
			if (c == '\r') {
				if (sb->sgetc() == '\n') sb->sbumpc();
				return true;
			}
			if (c == std::streambuf::traits_type::eof()) {
				if (line.empty()) {
					gin.setstate(std::ios::eofbit);
					return false;
				}
				return true;				
			}
			if (isprint(c)) line += (char)c;
		} while (true);
	}
};

struct oGZipFile {
	std::ofstream fout;
	boost::iostreams::filtering_ostream gout;

	oGZipFile(const std::string& output_file = "") {
		if (output_file != "") open(output_file);
	}

	oGZipFile(const oGZipFile& o) {
	}

	~oGZipFile() { close(); }

	bool open(const std::string& output_file) {
		size_t pos = output_file.find(".gz");
		bool isGZ = (pos != std::string::npos) && (pos + 3 == output_file.length());

		fout.open(output_file, std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
		if (isGZ) gout.push(boost::iostreams::gzip_compressor());
		gout.push(fout);

		return true;
	}

	bool close() {
		if (!gout.empty()) gout.reset();
		if (fout.is_open()) fout.close();
		return true;
	}

	void write(const Read& aread) {
		const std::string &outstr = aread.toString();
		gout.write(outstr.c_str(), outstr.length());
	}

	boost::iostreams::filtering_ostream& operator()() { return gout; }
};

#endif
