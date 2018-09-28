// Copyright (C) 2018 The Regents of the University of California (Regents).
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
// Author: Chris Sweeney (sweeney.chris.m@gmail.com)

#ifndef THEIA_MATH_RESERVOIR_SAMPLER_H_
#define THEIA_MATH_RESERVOIR_SAMPLER_H_

#include <vector>

#include <theia/util/random.h>

namespace theia {

// A reservoir sampler is a memory-constrained probabilistic sampler. It is
// suitable for drawing random samples from a sequence with a large or unknown
// number of items. As such, this class is suitable for drawing K random samples
// from a sequence of N >> K elements without having to hold all N elements in
// memory.
//
// The sample proceeds by choosing to keep the i-th element added to the sampler
// with a probability of K / i.
template <typename ElementType>
class ReservoirSampler {
 public:
  // TODO(sweeneychris): Add a constructor to set the RNG or the seed.

  // The number of elements we would like to sample from the entire sequence.
  explicit ReservoirSampler(const int num_elements_to_sample)
      : num_elements_to_sample_(num_elements_to_sample),
        num_elements_added_(0) {
    randomly_sampled_elements_.reserve(num_elements_to_sample_);
  }

  // Add a single element to the sampler. The element will be retained (i.e.
  // will be considered a part of the random K samples) with probability K / i,
  // where i is the # of elements added to the sampler so far.
  void AddElementToSampler(const ElementType& element) {
    // If we do not have enough samples yet, add the element to the sampling
    // with probabiliy of 1.
    if (num_elements_added_ < num_elements_to_sample_) {
      randomly_sampled_elements_.push_back(element);
    } else {
      // Otherwise, we want to add the new element to our sampling with if
      //   Rand(0.0, 1.0) < K / i.
      // This is equivalent to evaluating the probability that Random(0, N) < K
      // where N is the number of elements added so far, but this version avoids
      // costly division operators for each sample.
      const int modified_sample_probability =
          rng_.RandInt(0, num_elements_added_);
      if (modified_sample_probability < num_elements_to_sample_) {
        randomly_sampled_elements_[modified_sample_probability] = element;
      }
    }
    ++num_elements_added_;
  }

  // Returns all of the samples
  const std::vector<ElementType>& GetAllSamples() const {
    return randomly_sampled_elements_;
  }

  // Return the total number of elements added to the reservoir.
  int NumElementsAdded() const { return num_elements_added_; }

 private:
  // The number of elements we would like to sample from the entire sequence.
  const int num_elements_to_sample_;
  // The number of elements currently added to the sampler. This informs how to
  // probabilistically sample new data as it is added.
  int num_elements_added_;
  RandomNumberGenerator rng_;

  // The current random sampling of elements.
  std::vector<ElementType> randomly_sampled_elements_;
};
}  // namespace theia

#endif  // THEIA_MATH_RESERVOIR_SAMPLER_H_
