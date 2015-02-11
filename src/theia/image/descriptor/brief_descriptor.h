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

#ifndef THEIA_IMAGE_DESCRIPTOR_BRIEF_DESCRIPTOR_H_
#define THEIA_IMAGE_DESCRIPTOR_BRIEF_DESCRIPTOR_H_

#include <Eigen/Core>
#include <vector>

#include "theia/image/descriptor/binary_descriptor.h"
#include "theia/image/descriptor/descriptor_extractor.h"
#include "theia/image/keypoint_detector/keypoint.h"
#include "theia/util/util.h"

namespace theia {
template <class T> class Image;
typedef Image<float> FloatImage;

// A helper struct to keep track of the two pixel locations to sample in the
// BRIEF pattern.
struct BriefSamplePair {
  int pixel_points[4];
};

// BRIEF Descriptor as described by Michael Calonder, Vincent Lepetit, Christoph
// Strecha, and Pascal Fua, “BRIEF: Binary Robust Independent Elementary
// Features”, 11th European Conference on Computer Vision (ECCV), Heraklion,
// Crete. LNCS Springer, September 2010
//
// The idea behind the BRIEF descriptor is to random sample two points in a
// (blurred) patch. Bits in the binary descriptor are determined by comparing
// the intensities in the sample points.
class BriefDescriptorExtractor : public BinaryDescriptorExtractor {
 public:
  // Each BRIEF descriptor is extracted within a S x S patch centered around the
  // keypoint, where S is the patch_sample_size. The number of samples per
  // descriptor is equal to 8 * num_bytes.
  BriefDescriptorExtractor(const int patch_sample_size, const int num_bytes)
      : patch_sample_size_(patch_sample_size), num_samples_(8 * num_bytes) {}
  BriefDescriptorExtractor() : BriefDescriptorExtractor(48, 32) {}
  ~BriefDescriptorExtractor() {}

  // BRIEF must use the same pattern for all features that need to be compared
  // so the pattern is created with Initialize..
  bool Initialize();

  // Computes a descriptor at a single keypoint.
  bool ComputeDescriptor(const FloatImage& image,
                         const Keypoint& keypoint,
                         BinaryVectorX* descriptor);

  // Compute multiple descriptors for keypoints from a single image.
  bool ComputeDescriptors(
      const FloatImage& image,
      std::vector<Keypoint>* keypoints,
      std::vector<BinaryVectorX>* descriptors);

  bool DetectAndExtractDescriptors(
      const FloatImage& image,
      std::vector<Keypoint>* keypoints,
      std::vector<BinaryVectorX>* descriptors);

  const std::vector<BriefSamplePair>& SamplePairs() const;

 private:
  bool ComputeDescriptorInternal(const FloatImage& image,
                                 const Keypoint& keypoints,
                                 BinaryVectorX* descriptors);

  const int patch_sample_size_;
  const int num_samples_;
  std::vector<BriefSamplePair> pixel_samples_;

  DISALLOW_COPY_AND_ASSIGN(BriefDescriptorExtractor);
};

}  // namespace theia

#endif  // THEIA_IMAGE_DESCRIPTOR_BRIEF_DESCRIPTOR_H_
