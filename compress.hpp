#ifndef COMPRESS
#define COMPRESS

/** This header was inspired by fastp's (https://github.com/OpenGene/fastp) writer.h and writer.cpp */

#include <cstdio>
#include <cstdlib>
#include <libdeflate.h>

const size_t compressor_buffer_size = 1ULL << 12;
const double compression_ratio_upper_bound = 5.0; // 1:5


// Base class
struct Compressor {
	virtual ~Compressor() {}
	virtual bool needFlush(size_t size) const { return false; }
	virtual void writeToBuffer(char c) {}
	virtual void writeToBuffer(const char* data, size_t size) {}
	virtual size_t compress() { return 1; } // return number of bytes after compression
	virtual void flushOut(FILE* fo, size_t out_size) {} // write outBuffer to file
};


struct SingleThreadCompressor: public Compressor {
	size_t datsize, bufsize, out_bufsize;
	char *inBuffer;
	void *outBuffer;
	libdeflate_compressor* compressor;

	SingleThreadCompressor(size_t buffer_size = compressor_buffer_size, int compression_level = 6) {
		datsize = 0;
		bufsize = buffer_size;
		out_bufsize = size_t(buffer_size / compression_ratio_upper_bound + 0.5);

		inBuffer = (char*)malloc(buffer_size);
		if (inBuffer == NULL) {
			printf("Failed to allocate an inBuffer of size %zu in SingleThreadCompressor!\n", buffer_size);
			exit(-1);
		}

		compressor = libdeflate_alloc_compressor(compression_level);
		if (compressor == NULL) {
			printf("Failed to create a libdeflate compressor in SingleThreadCompressor!\n");
			exit(-1);
		}

		outBuffer = malloc(out_bufsize);
		if (outBuffer == NULL) {
			printf("Failed to allocate an outBuffer of size %zu in SingleThreadCompressor!\n", out_bufsize);
			exit(-1);			
		}
	}

	~SingleThreadCompressor() {
		libdeflate_free_compressor(compressor);
		free(inBuffer);
		free(outBuffer);
	}

	bool needFlush(size_t size) const {
		return datsize + size > bufsize;
	}

	void writeToBuffer(char c) {
		*(inBuffer + datsize) = c;
		++datsize;
	}

	void writeToBuffer(const char* data, size_t size) {
		memcpy(inBuffer + datsize, data, size);
		datsize += size;
	}

	size_t compress() {
		size_t bound, out_size;

		if (datsize == 0) return 0;

		bound = libdeflate_gzip_compress_bound(compressor, datsize);
		
		if (out_bufsize < bound) {
			outBuffer = realloc(outBuffer, bound);
			if (outBuffer == NULL) {
				printf("Failed to realloc outBuffer in compress function, SingleThreadCompressor!\n");
				exit(-1);
			}
			out_bufsize = bound;
		}

		out_size = libdeflate_gzip_compress(compressor, inBuffer, datsize, outBuffer, bound);
		if (out_size == 0) {
			printf("Failed to compress inBuffer in compress, SingleThreadCompressor!\n");
			exit(-1);
		}

		datsize = 0;

		return out_size;
	}

	void flushOut(FILE* fo, size_t out_size) {
		size_t written_size = fwrite(outBuffer, 1, out_size, fo);
		if (out_size != written_size) {
			printf("Cannot write the full compressed buffer to disk, flushOut, SingleThreadCompressor!\n");
			exit(-1);
		}
	}
};

#endif
