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
// Authors: Victor Fragoso (victor.fragoso@mail.wvu.edu)
//        : Benjamin Smith (bbsmith1@mix.wvu.edu)

#ifndef THEIA_MATCHING_GLOBAL_DESCRIPTOR_EXTRACTOR_H_
#define THEIA_MATCHING_GLOBAL_DESCRIPTOR_EXTRACTOR_H_

#include <vector>
#include <glog/logging.h>
#include <Eigen/Core>

namespace theia {
// This class computes a global descriptor for a given image.
//
// The use of an extractor should follow the next procedure:
//
// 1. After creation, initialize the instance by calling Initialize().
// 2. Train the extractor by calling Train().
// 3. Compute the global descriptor by calling Extract().
class GlobalDescriptorExtractor {
 public:
  GlobalDescriptorExtractor() {}
  virtual ~GlobalDescriptorExtractor() {}

  // Initializes the global descriptor extractor.
  virtual bool Initialize() {
    return true;
  }

  // Trains the instance. Returns true upon success and false otherwise.
  // Params:
  //   training_descriptors  The training local descriptors by using specific
  //     a TrainingDescriptorType (e.g., vector of Eigens (matrix)).
  virtual bool Train(const std::vector<Eigen::VectorXf>& training_descriptors) {
    LOG(FATAL) << "Invalid call to the descriptor extractor";
    return false;
  }

  // Computes the global descriptor from information about the image to
  // describe. The method returns true upon success and false otherwise.
  // Params:
  //   input_information  The input image information.
  //   global_descriptor  The computed global descriptor.
  virtual bool Extract(const std::vector<Eigen::VectorXf> input_information,
               Eigen::VectorXf* global_descriptor) const {
    LOG(FATAL) << "Invalid call to the descriptor extractor";
    return false;
  }

};

}  // namespace theia

#endif  // THEIA_MATCHING_GLOBAL_DESCRIPTOR_EXTRACTOR_H_
