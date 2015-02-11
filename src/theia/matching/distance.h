// Copyright (C) 2013 The Regents of the University of California (Regents).
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

#ifndef THEIA_MATCHING_DISTANCE_H_
#define THEIA_MATCHING_DISTANCE_H_

#include <Eigen/Core>
#include <glog/logging.h>

#include "theia/image/descriptor/binary_descriptor.h"

namespace theia {
// This file includes all of the distance metrics that are used:
// L2 distance for euclidean features.
// Hamming distance for binary vectors.
// TODO(cmsweeney): Could implement Hamming distance for FREAK (checks the first
// 128 bits, then checks the rest).

// Squared Euclidean distance functor. We let Eigen handle the SSE optimization.
// NOTE: This assumes that each vector has a unit norm:
//  ||x - y||^2 = ||x||^2 + ||y||^2 - 2*||x^t * y|| = 2 - 2 * x.dot(y).
struct L2 {
  typedef float DistanceType;
  typedef Eigen::VectorXf DescriptorType;

  DistanceType operator()(const Eigen::VectorXf& descriptor_a,
                          const Eigen::VectorXf& descriptor_b) const {
    DCHECK_EQ(descriptor_a.size(), descriptor_b.size());
    const DistanceType dist = 2.0 - 2.0 * descriptor_a.dot(descriptor_b);
    return dist;
  }
};

// Haming distance functor. We break the binary descriptor down into bytes and
// use a lookup table to make the xor really fast.
struct Hamming {
  typedef int DistanceType;
  typedef BinaryVectorX DescriptorType;

  DistanceType operator()(const BinaryVectorX& descriptor_a,
                          const BinaryVectorX& descriptor_b) const {
    DCHECK_EQ(descriptor_a.size(), descriptor_b.size());
    static constexpr unsigned char pop_count_table[] = {
      0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, 1, 2, 2, 3, 2, 3, 3, 4, 2,
      3, 3, 4, 3, 4, 4, 5, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3,
      3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3,
      4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 2, 3, 3, 4,
      3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5,
      6, 6, 7, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4,
      4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5,
      6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 2, 3, 3, 4, 3, 4, 4, 5,
      3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 3,
      4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 4, 5, 5, 6, 5, 6, 6, 7, 5, 6,
      6, 7, 6, 7, 7, 8
    };

    int result = 0;
    const unsigned char* char_a = descriptor_a.data();
    const unsigned char* char_b = descriptor_b.data();
    for (size_t i = 0;
         i < descriptor_a.size() / (8 * sizeof(unsigned char)); i++) {
      result += pop_count_table[char_a[i] ^ char_b[i]];
    }

    return result;
  }
};

}  // namespace theia

#endif  // THEIA_MATCHING_DISTANCE_H_
