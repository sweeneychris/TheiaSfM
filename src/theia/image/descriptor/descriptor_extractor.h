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

#ifndef THEIA_IMAGE_DESCRIPTOR_DESCRIPTOR_EXTRACTOR_H_
#define THEIA_IMAGE_DESCRIPTOR_DESCRIPTOR_EXTRACTOR_H_

#include <glog/logging.h>
#include <Eigen/Core>
#include <algorithm>
#include <vector>

#include "theia/image/keypoint_detector/keypoint.h"
#include "theia/util/util.h"

namespace theia {
class FloatImage;

// This interface class is meant to define all descriptor extractors. Different
// descriptor types may be easily implemented by deriving from this class.
class DescriptorExtractor {
 public:
  DescriptorExtractor() {}
  virtual ~DescriptorExtractor() {}

  // This method should be called before using any of the descriptor
  // extractors. The constuctor should only be used for get/set methods. Any
  // computationally expensive or non-trivial operations should go in the
  // Initialize method.
  virtual bool Initialize() { return true; }

  // Computes a floatdescriptor at a single keypoint.
  virtual bool ComputeDescriptor(const FloatImage& image,
                                 const Keypoint& keypoint,
                                 Eigen::VectorXf* descriptor) = 0;

  // Compute the descriptors for multiple keypoints in a given image. This
  // method will return all descriptors that could be extracted. If any
  // descriptors could not be extracted at a given keypoint, that keypoint will
  // be removed from the container. Returns true on success and false on
  // failure.
  virtual bool ComputeDescriptors(
      const FloatImage& image,
      std::vector<Keypoint>* keypoints,
      std::vector<Eigen::VectorXf>* descriptors);

  // Detects keypoints using the default method for the given descriptor. This
  // can be more efficient (e.g., with SIFT) because there is some overhead
  // required for creating the keypoint and descriptor objects.
  virtual bool DetectAndExtractDescriptors(
      const FloatImage& image,
      std::vector<Keypoint>* keypoints,
      std::vector<Eigen::VectorXf>* descriptors) = 0;

 private:
  DISALLOW_COPY_AND_ASSIGN(DescriptorExtractor);
};

}  // namespace theia

#endif  // THEIA_IMAGE_DESCRIPTOR_DESCRIPTOR_EXTRACTOR_H_
