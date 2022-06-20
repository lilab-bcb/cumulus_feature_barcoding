/**
 * Slightly modified by Bo Li from Rob's FQfeeder: https://github.com/rob-p/FQFeeder. See below for the original LICENSE.
 * 
**/

/**
Copyright (c) 2016, Rob Patro
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

* Neither the name of FQFeeder nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
**/

#ifndef THREAD_UTILS_HPP
#define THREAD_UTILS_HPP

#include <cassert>
#include <chrono>
#include <pthread.h>
#include <random>
#include <thread>

#if defined(__SSE2__)
 #if defined(HAVE_SIMDE)
  #include "simde/x86/sse2.h"
 #else
  #include <emmintrin.h>
 #endif
#endif

// Most of this code is taken directly from
// https://github.com/geidav/spinlocks-bench/blob/master/os.hpp. However, things
// may be renamed, modified, or randomly mangled over time.
#define ALWAYS_INLINE inline __attribute__((__always_inline__))

static const constexpr size_t MIN_BACKOFF_ITERS = 32;
static const size_t MAX_BACKOFF_ITERS = 1024;

ALWAYS_INLINE static void cpuRelax() {
#if defined(__SSE2__)  // AMD and Intel
  #if defined(HAVE_SIMDE)
    simde_mm_pause();
  #else
    _mm_pause();
  #endif
#elif defined(__i386__) || defined(__x86_64__)
  asm volatile("pause");
#elif defined(__aarch64__)
  asm volatile("wfe");
#elif defined(__armel__) || defined(__ARMEL__)
  asm volatile ("nop" ::: "memory");
#elif defined(__arm__) || defined(__aarch64__)
  __asm__ __volatile__ ("yield" ::: "memory");
#elif defined(__ia64__)  // IA64
  __asm__ __volatile__ ("hint @pause");
#elif defined(__powerpc__) || defined(__ppc__) || defined(__PPC__)
   __asm__ __volatile__ ("or 27,27,27" ::: "memory");
#else  // everything else.
   asm volatile ("nop" ::: "memory");
#endif
}

ALWAYS_INLINE void yieldSleep() {
  using namespace std::chrono;
  std::chrono::microseconds ytime(500);
  std::this_thread::sleep_for(ytime);
}

ALWAYS_INLINE void backoffExp(size_t& curMaxIters) {
  thread_local std::uniform_int_distribution<size_t> dist;

  // see : https://github.com/coryan/google-cloud-cpp-common/blob/a6e7b6b362d72451d6dc1fec5bc7643693dbea96/google/cloud/internal/random.cc
  #if defined(__linux) && defined(__GLIBCXX__) && __GLIBCXX__ >= 20200128
    thread_local std::random_device rd("/dev/urandom");
  #else
    thread_local std::random_device rd;
  #endif  // defined(__GLIBCXX__) && __GLIBCXX__ >= 20200128

  thread_local std::minstd_rand gen(rd());
  const size_t spinIters =
      dist(gen, decltype(dist)::param_type{0, curMaxIters});
  curMaxIters = std::min(2 * curMaxIters, MAX_BACKOFF_ITERS);
  for (size_t i = 0; i < spinIters; i++) {
    cpuRelax();
  }
}

ALWAYS_INLINE void backoffOrYield(size_t& curMaxDelay) {
  if (curMaxDelay >= MAX_BACKOFF_ITERS) {
    yieldSleep();
    curMaxDelay = MIN_BACKOFF_ITERS;
  }
  backoffExp(curMaxDelay);
}

#endif 
