#ifndef GZIP_UTILS
#define GZIP_UTILS

#include <cctype>
#include <string>
#include <istream>
#include <streambuf>
#include <fstream>

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
    
	iGZipKseqFile(const std::string& input_file) {
		sequence_file = gzopen(input_file.c_str(), "r");
		if (sequence_file == NULL) {
			printf("Cannot find sequence file %s!\n", input_file.c_str());
			exit(-1);
		}
		sequence_kseq = kseq_init(sequence_file);
	}

	~iGZipKseqFile() {
		kseq_destroy(sequence_kseq);
		gzclose(sequence_file);
		sequence_kseq = NULL;
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

#endif
