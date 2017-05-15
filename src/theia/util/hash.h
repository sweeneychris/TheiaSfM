// Copyright (C) 2014 The Regents of the University of California (Regents).
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//
//     * Neither the name of The Regents or University of California nor the
//       names of its contributors may be used to endorse or promote products
//       derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Please contact the author of this library if you have any questions.
// Author: Chris Sweeney (cmsweeney@cs.ucsb.edu)

#ifndef THEIA_UTIL_HASH_H_
#define THEIA_UTIL_HASH_H_

#include <Eigen/Core>
#include <utility>

// This file defines hash functions for stl containers.
namespace std {
namespace {

// Combines the hash of v with the current hash value seed. This is the
// recommended approach from Boost.
template <class T>
inline void HashCombine(const T& v, std::size_t* seed) {
  std::hash<T> hasher;
  *seed ^= hasher(v) + 0x9e3779b9 + (*seed << 6) + (*seed >> 2);
}

}  // namespace

// STL does not implement hashing for pairs, so a simple pair hash is done here.
template <typename T1, typename T2> struct hash<std::pair<T1, T2> > {
 public:
  size_t operator()(const std::pair<T1, T2>& e) const {
    size_t seed = 0;
    HashCombine(e.first, &seed);
    HashCombine(e.second, &seed);
    return seed;
  }
};

// A generic hash function that hashes constant size Eigen matrices and vectors.
template <typename T, int N, int M>
struct hash<Eigen::Matrix<T, N, M> > {
  size_t operator()(const Eigen::Matrix<T, N, M>& matrix) const {
    size_t seed = 0;
    const T* data = matrix.data();
    for (size_t i = 0; i < matrix.size(); ++i) {
      seed ^= std::hash<T>()(data[i]) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    }
    return seed;
  }
};

}  // namespace std

#endif  // THEIA_UTIL_HASH_H_
