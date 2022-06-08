#ifndef GZIP_UTILS
#define GZIP_UTILS

#include <cctype>
#include <cassert>
#include <string>
#include <cstdio>

#include "external/slw287r_trimadap/izlib.h"
#include "external/kseq.h"

struct Read {
	std::string name, comment, seq, qual; // qual.empty() -> FASTQ reads

	std::string toString() const {
		return "@" + name + (comment.empty() ? "" : (" " + comment)) + "\n" + seq + "\n+\n" + qual + "\n";
	}
};

struct iGZipFile {
	KSEQ_INIT(gzFile, gzread);

	gzFile sequence_file = NULL;
	kseq_t* sequence_kseq = NULL;


	iGZipFile(const std::string& input_file) {
		sequence_file = gzopen(input_file.c_str(), "r");
		if (sequence_file == NULL) {
			printf("Cannot find sequence file %s!\n", input_file.c_str());
			exit(-1);
		}
		sequence_kseq = kseq_init(sequence_file);
	}

	~iGZipFile() {
		kseq_destroy(sequence_kseq);
		gzclose(sequence_file);
		sequence_kseq = NULL;
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
			printf("Error reading stream; didn't reach the end of sequence file, which might be corrupted!\n");
			exit(-1);			
		}
		else if (length == -2) {
			printf("Truncated quality string!\n");
			exit(-1);			
		}
		else {
			assert(length == 0);
			printf("Detected a read with 0 sequence length!\n");
			exit(-1);
		}

		return true;
	}


	bool next(std::string& line) { // Only get one line
		int ret = ks_getuntil(sequence_kseq->f, KS_SEP_LINE, &sequence_kseq->comment, 0);

		if (ret == -3) {
			printf("Error reading stream in next(line)!\n");
			exit(-1);						
		}

		if (ret >= 0) {
			line = sequence_kseq->comment.s;
			return true;
		}

		return false;
	}
};

#endif
