//
// Created by Ivan Radkevich on 6/13/22.
// Rewrote by Bo Li on 6/18/22.
// This is a simplified reimplementation of Rob Patro's FQFeeder: https://github.com/rob-p/FQFeeder [BSD-3 license]
//

#ifndef READPARSER_HPP
#define READPARSER_HPP

#include <atomic>
#include <vector>
#include <cstdio>
#include <cassert>
#include <memory>

#include "gzip_utils.hpp"
#include "external/concurrentqueue.h"
#include "external/thread_utils.hpp"


struct ReadTuple {
    std::vector<Read> mates;

    ReadTuple(int n_mates) : mates(n_mates) {}
    Read& operator[](int i) { return mates[i]; }
};

struct ReadChunk {
    size_t nreads; // how many reads we currently have?
    std::vector<ReadTuple> reads;

    ReadChunk(size_t size, int n_mates) : reads(size, ReadTuple(n_mates)) {}
    ReadTuple& operator[](size_t i) { return reads[i]; }
    std::vector<ReadTuple>::iterator begin() { return reads.begin(); }
    std::vector<ReadTuple>::iterator end() { return reads.begin() + nreads; }
};

struct ReadGroup {
    std::unique_ptr<ReadChunk> chunk_{nullptr};
    moodycamel::ProducerToken ptSpace;
    moodycamel::ConsumerToken ctRead;

    ReadGroup(moodycamel::ProducerToken&& pt, moodycamel::ConsumerToken&& ct) : ptSpace(std::move(pt)), ctRead(std::move(ct)) {}

    bool empty() const { return chunk_.get() == nullptr; }

    ReadTuple& operator[](size_t i) { return (*chunk_)[i]; };
    std::vector<ReadTuple>::iterator begin() { return chunk_->begin(); }
    std::vector<ReadTuple>::iterator end() { return chunk_->end(); }
};

class ReadParser {
public:
    ReadParser(std::vector<std::vector<std::string>>& input_files, int numConsumers, int numParsers = 1, size_t chunkSize = 100000) : chunkSize_(chunkSize), input_files_(input_files) {
        int n_files = input_files.size();
        n_mates_ = input_files[0].size();

        if (numParsers > n_files) {
            printf("Detected more parsers than number of files; setting numParsers to %d.\n", n_files);
            numParsers = n_files;
        }
        
        // Initialize queues
        fileQueue_ = moodycamel::ConcurrentQueue<int>(numParsers); // Only need one producer
        moodycamel::ProducerToken ptFile(fileQueue_);
        for (int i = 0; i < n_files; ++i) assert(fileQueue_.enqueue(ptFile, i));

        spaceQueue_ = moodycamel::ConcurrentQueue<std::unique_ptr<ReadChunk>>(4 * numConsumers, 1 + numConsumers, 0); // blocks of empty space (ReadChunk) to fill in, the extra one is for this thread
        readQueue_ = moodycamel::ConcurrentQueue<std::unique_ptr<ReadChunk>>(4 * numConsumers, numParsers, 0); // blocks of filled in reads

        // Each parser should have a consumer token for spaceQueue_ to get empty spaces and one producer token for the readQueue_ to load reads
        for (int i = 0; i < numParsers; ++i) {
            ctSpaceQueue_.emplace_back(new moodycamel::ConsumerToken(spaceQueue_));
            ptReadQueue_.emplace_back(new moodycamel::ProducerToken(readQueue_));
        }

        // create empty chunks and push to spaceQueue_
        moodycamel::ProducerToken ptoken(spaceQueue_);
        for (int i = 0; i < 4 * numConsumers; ++i) {
          auto chunk = std::unique_ptr<ReadChunk>(new ReadChunk(chunkSize_, n_mates_));
          assert(spaceQueue_.enqueue(ptoken, std::move(chunk)));
        }

        //start parsing
        numParsing_ = 0;
        for (int i = 0; i < numParsers; ++i) {
          ++numParsing_;
          parsingThreads_.emplace_back(new std::thread([this, i]() {
                this->parse_read_tuples(this->ctSpaceQueue_[i].get(), this->ptReadQueue_[i].get());
            }));
        }
    }

    ~ReadParser() {
        for (auto& thread: parsingThreads_) thread->join();
    }

    ReadGroup getReadGroup() {
        return ReadGroup(moodycamel::ProducerToken(spaceQueue_), moodycamel::ConsumerToken(readQueue_));
    }

    bool refill(ReadGroup& rg) {
        if (!rg.empty()) {
            assert(spaceQueue_.enqueue(rg.ptSpace, std::move(rg.chunk_)));
            assert(rg.empty());
        }
        
        size_t curMaxDelay = MIN_BACKOFF_ITERS;
        while (numParsing_ > 0) {
            if (readQueue_.try_dequeue(rg.ctRead, rg.chunk_)) return true;
            backoffOrYield(curMaxDelay);
        }

        return readQueue_.try_dequeue(rg.ctRead, rg.chunk_);
    }

private:
    size_t chunkSize_;
    int n_mates_;
    std::vector<std::vector<std::string>> input_files_;

    std::atomic<int> numParsing_;

    moodycamel::ConcurrentQueue<int> fileQueue_; // hold file positions in input_files_
    moodycamel::ConcurrentQueue<std::unique_ptr<ReadChunk>> spaceQueue_, readQueue_;

    std::vector<std::unique_ptr<moodycamel::ProducerToken>> ptReadQueue_;
    std::vector<std::unique_ptr<moodycamel::ConsumerToken>> ctSpaceQueue_;

    std::vector<std::unique_ptr<std::thread>> parsingThreads_;


    bool load_one_tuple(std::vector<iGZipFile>& input_streams, ReadTuple *rt) {
        int cnt = 0;
        for (int i = 0; i < n_mates_; ++i)
            cnt += input_streams[i].next((*rt)[i]);
        if (cnt > 0 && cnt != n_mates_) {
            printf("Detected mate files with different number of lines!\n");
            exit(-1);
        }
        return cnt > 0;
    }

    void parse_read_tuples(moodycamel::ConsumerToken* ctSpace, moodycamel::ProducerToken* ptRead) {
        int fid;
        size_t cur_size, curMaxDelay;
        ReadTuple *rt;
        std::unique_ptr<ReadChunk> buffer;
        std::vector<iGZipFile> input_streams;
        
        // load reading buffer
        curMaxDelay = MIN_BACKOFF_ITERS;
        while (!spaceQueue_.try_dequeue(*ctSpace, buffer)) backoffOrYield(curMaxDelay);

        cur_size = 0;
        rt = &((*buffer)[cur_size]);

        while (fileQueue_.try_dequeue(fid)) {
            // prepare input gzip "streams"
            input_streams.clear();
            for (int i = 0; i < n_mates_; ++i) input_streams.emplace_back(input_files_[fid][i]);

            // load read tuples
            while (load_one_tuple(input_streams, rt)) {
                ++cur_size;
                if (cur_size == chunkSize_) {
                    buffer->nreads = cur_size;
                    curMaxDelay = MIN_BACKOFF_ITERS;
                    while (!readQueue_.try_enqueue(*ptRead, std::move(buffer))) backoffOrYield(curMaxDelay);
                    cur_size = 0;
                    curMaxDelay = MIN_BACKOFF_ITERS;
                    while (!spaceQueue_.try_dequeue(*ctSpace, buffer)) backoffOrYield(curMaxDelay);
                }
                rt = &((*buffer)[cur_size]);
            }
        }

        if (cur_size > 0) {
            buffer->nreads = cur_size;
            curMaxDelay = MIN_BACKOFF_ITERS;
            while (!readQueue_.try_enqueue(*ptRead, std::move(buffer))) backoffOrYield(curMaxDelay);
        }
        
        --numParsing_;
    }
};

#endif
