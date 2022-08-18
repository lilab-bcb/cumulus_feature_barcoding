#ifndef COMPRESS
#define COMPRESS

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <thread>
#include <libdeflate.h>

const size_t compressor_buffer_size = 1ULL << 23; // 8M buffer

const size_t bgzf_block_size = 0xff00; // same as htslib
const size_t bgzf_block_header_length = 18;
const size_t bgzf_block_footer_length = 8;

// The annotation below is from htslib
/* BGZF/GZIP header (specialized from RFC 1952; little endian):
 +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
 | 31|139|  8|  4|              0|  0|255|      6| 66| 67|      2|BLK_LEN|
 +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
  BGZF extension:
                ^                              ^   ^   ^
                |                              |   |   |
               FLG.EXTRA                     XLEN  B   C

  BGZF format is compatible with GZIP. It limits the size of each compressed
  block to 2^16 bytes and adds and an extra "BC" field in the gzip header which
  records the size.

*/
const uint8_t bgzf_block_header[19] = "\037\213\010\4\0\0\0\0\0\377\6\0\102\103\2\0\0\0";
const uint8_t bgzf_empty_block[29] = "\037\213\010\4\0\0\0\0\0\377\6\0\102\103\2\0\033\0\3\0\0\0\0\0\0\0\0\0";


// Base class
struct Compressor {
	virtual ~Compressor() {}
	virtual bool needFlush(size_t size) const { return false; }
	virtual void writeToBuffer(char c) {}
	virtual void writeToBuffer(const char* data, size_t size) {}
	virtual size_t compress() { return 1; } // return number of bytes after compression
	virtual void flushOut(FILE* fo, size_t out_size, bool last_flush = false) {} // write outBuffer to file
};


struct SingleThreadCompressor: public Compressor {
	size_t datsize, bufsize, out_bufsize;
	int compression_level;
	uint8_t *inBuffer, *outBuffer;
	libdeflate_compressor* compressor;
	bool is_bgzf;

	SingleThreadCompressor(size_t buffer_size = compressor_buffer_size, int compression_level = 6, bool is_bgzf = false) : datsize(0), bufsize(buffer_size), compression_level(compression_level), is_bgzf(is_bgzf) {
		inBuffer = (uint8_t*)malloc(buffer_size);
		if (inBuffer == NULL) {
			printf("Failed to allocate an inBuffer of size %zu in SingleThreadCompressor!\n", buffer_size);
			exit(-1);
		}

		compressor = libdeflate_alloc_compressor(compression_level);
		if (compressor == NULL) {
			printf("Failed to create a libdeflate compressor in SingleThreadCompressor!\n");
			exit(-1);
		}

		out_bufsize = 0;
		if (!is_bgzf) {
			out_bufsize = libdeflate_gzip_compress_bound(compressor, buffer_size);
		} 
		else {
			size_t block_size = libdeflate_deflate_compress_bound(compressor, bgzf_block_size) + bgzf_block_header_length + bgzf_block_footer_length;
			out_bufsize = (buffer_size / bgzf_block_size) * block_size;
			size_t rest_size = buffer_size % bgzf_block_size;
			if (rest_size > 0) {
				out_bufsize += libdeflate_deflate_compress_bound(compressor, rest_size) + bgzf_block_header_length + bgzf_block_footer_length;
			}
		}

		outBuffer = (uint8_t*)malloc(out_bufsize);
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
		size_t out_size = 0;
		
		if (datsize == 0) return out_size;

		if (!is_bgzf) {
			out_size = libdeflate_gzip_compress(compressor, inBuffer, datsize, outBuffer, out_bufsize);
			if (out_size == 0) {
				printf("Failed to compress inBuffer in compress, SingleThreadCompressor!\n");
				exit(-1);
			}
		}
		else {
			size_t in_nbytes = 0, out_block_size = 0, bsize;
			uint8_t *in_buff = nullptr, *out_buff = nullptr;
			uint32_t crc;

			for (size_t i = 0; i < datsize; i += bgzf_block_size) {
				in_buff = inBuffer + i;
				out_buff = outBuffer + out_size;

				in_nbytes = datsize - i;
				if (in_nbytes > bgzf_block_size) in_nbytes = bgzf_block_size;

				// write header
				memcpy(out_buff, bgzf_block_header, bgzf_block_header_length);
				
				if (compression_level == 0) {
					out_buff[bgzf_block_header_length] = 1; // BFINAL=1, BTYPE=00; see RFC1951; from htslib
					out_buff[bgzf_block_header_length + 1] = in_nbytes & 0xff;
					out_buff[bgzf_block_header_length + 2] = (in_nbytes >> 8) & 0xff;
					out_buff[bgzf_block_header_length + 3] = (~in_nbytes) & 0xff;
					out_buff[bgzf_block_header_length + 4] = ((~in_nbytes) >> 8) & 0xff;
					memcpy(out_buff + bgzf_block_header_length + 5, in_buff, in_nbytes);
					out_block_size = in_nbytes + 5 + bgzf_block_header_length + bgzf_block_footer_length;
				}
				else {
					out_block_size = libdeflate_deflate_compress(compressor, in_buff, in_nbytes, out_buff + bgzf_block_header_length, out_bufsize - out_size - bgzf_block_header_length);
					if (out_block_size == 0) {
						printf("Failed to compress BGZF block %d in compress, SingleThreadCompressor!\n", int(i / bgzf_block_size));
						exit(-1);
					}
					out_block_size += bgzf_block_header_length + bgzf_block_footer_length;
				}

				// fill in BSIZE
				bsize = out_block_size - 1;
				out_buff[bgzf_block_header_length - 2] = bsize & 0xff;
				out_buff[bgzf_block_header_length - 1] = (bsize >> 8) & 0xff;

				// writer footer
				out_buff += out_block_size - bgzf_block_footer_length;

				crc = libdeflate_crc32(0, in_buff, in_nbytes);
				
				out_buff[0] = crc & 0xff;
				out_buff[1] = (crc >> 8) & 0xff;
				out_buff[2] = (crc >> 16) & 0xff;
				out_buff[3] = (crc >> 24) & 0xff;				

				out_buff[4] = in_nbytes & 0xff;
				out_buff[5] = (in_nbytes >> 8) & 0xff;
				out_buff[6] = (in_nbytes >> 16) & 0xff;
				out_buff[7] = (in_nbytes >> 24) & 0xff;

				out_size += out_block_size;
			}
		}

		datsize = 0;

		return out_size;
	}

	void flushOut(FILE* fo, size_t out_size, bool last_flush = false) {
		size_t written_size;

		if (out_size > 0) {
			written_size = fwrite(outBuffer, 1, out_size, fo);
			if (written_size != out_size) {
				printf("Cannot write the full compressed buffer to disk, flushOut, SingleThreadCompressor!\n");
				exit(-1);
			}
			out_size = 0;
		}

		if (last_flush && is_bgzf) { // write last empty block for BGZF
			memcpy(outBuffer, bgzf_empty_block, 28);
			written_size = fwrite(outBuffer, 1, 28, fo);
			if (written_size != 28) {
				printf("Cannot write the empty BGZF block to disk, flushOut, SingleThreadCompressor!\n");
				exit(-1);				
			}
		}
	}
};


void perform_compression(SingleThreadCompressor* compressor, size_t* out_size) {
	*out_size = compressor->compress();
}

struct MultiThreadsCompressor: public Compressor {
	int num_threads, cur_pos; // current compressor to feed data
	SingleThreadCompressor **compressors;
	size_t *out_sizes;
	std::thread **threads;
	bool is_bgzf;

	MultiThreadsCompressor(int num_threads, size_t buffer_size = compressor_buffer_size, int compression_level = 6, bool is_bgzf = false) : num_threads(num_threads), cur_pos(0), is_bgzf(is_bgzf) {
		compressors = new SingleThreadCompressor*[num_threads];
		out_sizes = new size_t[num_threads]();
		threads = new std::thread*[num_threads];
		for (int i = 0; i < num_threads; ++i) {
			compressors[i] = new SingleThreadCompressor(buffer_size, compression_level, is_bgzf);
		}
	}

	~MultiThreadsCompressor() {
		for (int i = 0; i < num_threads; ++i) {
			delete compressors[i];
		}
		delete[] compressors;
		delete[] out_sizes;
		delete[] threads;
	}

	bool needFlush(size_t size) const {
		return (cur_pos + 1 == num_threads) && compressors[cur_pos]->needFlush(size);
	}

	void writeToBuffer(char c) {
		if (compressors[cur_pos]->needFlush(1)) ++cur_pos;
		compressors[cur_pos]->writeToBuffer(c);
	}

	void writeToBuffer(const char* data, size_t size) {
		if (compressors[cur_pos]->needFlush(size)) ++cur_pos;
		compressors[cur_pos]->writeToBuffer(data, size);
	}

	size_t compress() {
		int n_active = cur_pos + (compressors[cur_pos]->datsize > 0);
		cur_pos = 0;

		if (n_active == 0) return 0;
		if (n_active == 1) {
			out_sizes[0] = compressors[0]->compress();
			return out_sizes[0];
		}

		for (int i = 0; i < n_active; ++i)
			threads[i] = new std::thread(perform_compression, compressors[i], &out_sizes[i]);
		for (int i = 0; i < n_active; ++i) threads[i]->join();
		for (int i = 0; i < n_active; ++i) delete threads[i];

		size_t out_size = 0;
		for (int i = 0; i < n_active; ++i) out_size += out_sizes[i];

		return out_size;
	}

	void flushOut(FILE* fo, size_t out_size, bool last_flush = false) {
		int i = 0;
		while (i < num_threads && out_sizes[i] > 0) {
			compressors[i]->flushOut(fo, out_sizes[i]);
			out_sizes[i] = 0;
			++i;
		}

		if (last_flush && is_bgzf)
			compressors[0]->flushOut(fo, 0, true);
	}
};

#endif
