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

#ifndef THEIA_IMAGE_DESCRIPTOR_FREAK_DESCRIPTOR_H_
#define THEIA_IMAGE_DESCRIPTOR_FREAK_DESCRIPTOR_H_

#include <vector>

#include "theia/image/descriptor/binary_descriptor.h"
#include "theia/image/descriptor/descriptor_extractor.h"
#include "theia/image/keypoint_detector/keypoint.h"
#include "theia/util/util.h"

namespace theia {
template <class T> class Image;
typedef Image<float> FloatImage;

class FreakDescriptorExtractor : public BinaryDescriptorExtractor {
 public:
  // Params:
  //  rotation_invariant: Enable orientation normalization.
  //  scale_invariant: Enable scale normalization.
  //  num_octaves: Number of octaves covered by the keypoints.
  FreakDescriptorExtractor(bool rotation_invariant,
                           bool scale_invariant,
                           int num_octaves)
      : rotation_invariant_(rotation_invariant),
        scale_invariant_(scale_invariant),
        num_octaves_(num_octaves),
        pattern_scale_(22.0) {}

  FreakDescriptorExtractor() : FreakDescriptorExtractor(true, true, 4) {}
  ~FreakDescriptorExtractor() {}

  // Initializes the sampling patterns and local variables.
  bool Initialize();

  // Computes a descriptor at a single keypoint.
  bool ComputeDescriptor(const FloatImage& image,
                         const Keypoint& keypoint,
                         BinaryVectorX* descriptor);

  // Compute multiple descriptors for keypoints from a single image. Note this
  // may return null for some of the descriptors if they cannot be computed!
  // Typically this only happens when it is too close to the border.
  bool ComputeDescriptors(const FloatImage& image,
                          std::vector<Keypoint>* keypoints,
                          std::vector<BinaryVectorX>* descriptors);

  bool DetectAndExtractDescriptors(
    const FloatImage& image,
    std::vector<Keypoint>* keypoints,
    std::vector<BinaryVectorX>* descriptors);

 private:
  uchar MeanIntensity(const Image<uchar>& image, const Image<int>& integral,
                      const float kp_x, const float kp_y,
                      const unsigned int scale, const unsigned int rot,
                      const unsigned int point) const;

  static const int kNumScales_ = 64;
  static const int kNumPairs_ = 512;
  static const int kNumOrientationPairs_ = 45;

  bool rotation_invariant_;
  bool scale_invariant_;
  bool num_octaves_;
  const float pattern_scale_;

  std::vector<int> selected_pairs_;

  struct PatternPoint {
    // x coordinate relative to center
    float x;
    // x coordinate relative to center
    float y;
    // Gaussian smoothing sigma
    float sigma;
  };

  struct DescriptionPair {
    // index of the first point
    unsigned char i;
    // index of the second point
    unsigned char j;
  };

  struct OrientationPair {
    // index of the first point
    unsigned char i;
    // index of the second point
    unsigned char j;
    // dx/(norm_sq))*4096
    int weight_dx;
    // dy/(norm_sq))*4096
    int weight_dy;
  };

  // look-up table for the pattern points (position+sigma of all points at all
  // scales and orientation)
  std::vector<PatternPoint> pattern_lookup_;
  // size of the pattern at a specific scale (used to check if a point is within
  // image boundaries)
  int pattern_sizes_[kNumScales_];
  DescriptionPair description_pairs_[kNumPairs_];
  OrientationPair orientation_pairs_[kNumOrientationPairs_];

  DISALLOW_COPY_AND_ASSIGN(FreakDescriptorExtractor);
};

}  // namespace theia

#endif  // THEIA_IMAGE_DESCRIPTOR_FREAK_DESCRIPTOR_H_
