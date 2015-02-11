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

#ifndef THEIA_IMAGE_DESCRIPTOR_BRISK_DESCRIPTOR_H_
#define THEIA_IMAGE_DESCRIPTOR_BRISK_DESCRIPTOR_H_

#include <Eigen/Core>
#include <vector>

#include "theia/image/descriptor/binary_descriptor.h"
#include "theia/image/descriptor/descriptor_extractor.h"
#include "theia/image/keypoint_detector/brisk_detector.h"
#include "theia/image/keypoint_detector/keypoint.h"
#include "theia/util/util.h"

namespace theia {
template <class T> class Image;
typedef Image<float> FloatImage;

// BRISK Descriptor ported from reference code of Stefan Leutenegger, Margarita
// Chli and Roland Siegwart, "BRISK: Binary Robust Invariant Scalable
// Keypoints", in Proceedings of the IEEE International Conference on Computer
// Vision (ICCV 2011). NOTE: because this code is ported, it does not adhere to
// the same style guides as the rest of the code.
class BriskDescriptorExtractor : public BinaryDescriptorExtractor {
 public:
  // Set rotation_invariant and scale_invariant to true to calculate descriptors
  // that are invariant to under those shifts.
  BriskDescriptorExtractor(bool rotation_invariant,
                           bool scale_invariant,
                           float pattern_scale);
  BriskDescriptorExtractor() : BriskDescriptorExtractor(true, true, 1.0) {}
  ~BriskDescriptorExtractor();

  // Computes a descriptor at a single keypoint.
  bool ComputeDescriptor(const FloatImage& image,
                         const Keypoint& keypoint,
                         BinaryVectorX* descriptor);

  // Compute multiple descriptors for keypoints from a single image.
  bool ComputeDescriptors(
      const FloatImage& image,
      std::vector<Keypoint>* keypoints,
      std::vector<BinaryVectorX>* descriptors);

  // Detect keypoints using the Sift keypoint detector and extracts them at the
  // same time.
  bool DetectAndExtractDescriptors(
      const FloatImage& image,
      std::vector<Keypoint>* keypoints,
      std::vector<BinaryVectorX>* descriptors);

 private:
  // Call this to generate the kernel:
  // circle of radius r (pixels), with n points;
  // short pairings with dMax, long pairings with dMin
  void generateKernel(std::vector<float>* radiusList,
                      std::vector<int>* numberList, float dMax = 5.85f,
                      float dMin = 8.2f,
                      std::vector<int> indexChange = std::vector<int>());

  inline int smoothedIntensity(const Image<unsigned char>& image,
                               const Image<int>& integral, const float key_x,
                               const float key_y, const unsigned int scale,
                               const unsigned int rot,
                               const unsigned int point) const;

  bool rotation_invariance_;
  bool scale_invariance_;

  // Pattern Properties.
  // some helper structures for the Brisk pattern representation
  struct BriskPatternPoint {
    float x;      // x coordinate relative to center
    float y;      // x coordinate relative to center
    float sigma;  // Gaussian smoothing sigma
  };
  struct BriskShortPair {
    unsigned int i;  // index of the first pattern point
    unsigned int j;  // index of other pattern point
  };
  struct BriskLongPair {
    unsigned int i;   // index of the first pattern point
    unsigned int j;   // index of other pattern point
    int weighted_dx;  // 1024.0/dx
    int weighted_dy;  // 1024.0/dy
  };

  BriskPatternPoint* pattern_points_;
  unsigned int points_;
  // lists the scaling per scale index [scale]
  float* scale_list_;
  // lists the total pattern size per scale index [scale]
  unsigned int* size_list_;
  // scales discretization
  static const unsigned int scales_;
  // span of sizes 40->4 Octaves - else, this needs to be adjusted...
  static const float scale_range_;
  // discretization of the rotation look-up
  static const unsigned int n_rot_;

  // pairs
  // number of uchars the descriptor consists of
  int strings_;
  // short pair maximum distance
  float dMax_;
  // long pair maximum distance
  float dMin_;
  // d<_dMax
  BriskShortPair* short_pairs_;
  // d>_dMin
  BriskLongPair* long_pairs_;
  // number of shortParis
  unsigned int no_short_pairs_;
  // number of longParis
  unsigned int no_long_pairs_;

  // general
  static const float basic_size_;

  DISALLOW_COPY_AND_ASSIGN(BriskDescriptorExtractor);
};

}  // namespace theia

#endif  // THEIA_IMAGE_DESCRIPTOR_BRISK_DESCRIPTOR_H_
