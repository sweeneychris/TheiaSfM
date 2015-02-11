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

#include "theia/image/descriptor/brief_descriptor.h"

#include <Eigen/Core>
#include <algorithm>
#include <bitset>
#include <vector>

#include "theia/image/descriptor/binary_descriptor.h"
#include "theia/image/image.h"
#include "theia/image/descriptor/descriptor_extractor.h"
#include "theia/image/keypoint_detector/brisk_detector.h"
#include "theia/image/keypoint_detector/keypoint.h"
#include "theia/util/random.h"
#include "theia/util/util.h"

namespace theia {

namespace {

inline int ClampValue(const int min, const int max, const int val) {
  return std::max(std::min(val, max), min);
}

bool KeypointIsTooCloseToBorder(const Keypoint& keypoint,
                                const int image_rows,
                                const int image_cols,
                                const int patch_sample_size) {
  const int min_pixel_point = -patch_sample_size / 2;
  const int max_pixel_point = min_pixel_point + patch_sample_size - 1;

  return keypoint.x() + min_pixel_point < 0 ||
         keypoint.x() + max_pixel_point >= image_cols ||
         keypoint.y() + min_pixel_point < 0 ||
         keypoint.y() + max_pixel_point >= image_rows;
}
}  // namespace

// Generate pixels to sample.
bool BriefDescriptorExtractor::Initialize() {
  CHECK_GT(num_samples_, 0);
  CHECK_GT(patch_sample_size_, 0);
  CHECK_EQ(num_samples_ % 128, 0)
      << "BRIEF currently only supports num_samples being a multiple of 128";

  // The boundaries of the patch to sample.
  const int min_pixel_point = -patch_sample_size_ / 2;
  const int max_pixel_point = min_pixel_point + patch_sample_size_ - 1;

  // Using a zero-mean guassian distribution with a sigma of S * S / 25 is
  // reported to give good results in the paper.
  const double sample_sigma = patch_sample_size_ * patch_sample_size_ / 25.0;

  InitRandomGenerator();
  pixel_samples_.reserve(num_samples_);
  for (int i = 0; i < num_samples_; i++) {
    BriefSamplePair brief_sample_pair;
    for (int j = 0; j < 4; j++) {
      brief_sample_pair.pixel_points[j] =
          ClampValue(min_pixel_point, max_pixel_point,
                     static_cast<int>(RandGaussian(0.0, sample_sigma)));
    }
    pixel_samples_.emplace_back(brief_sample_pair);
  }

  return true;
}

bool BriefDescriptorExtractor::ComputeDescriptor(
    const FloatImage& image,
    const Keypoint& keypoint,
    BinaryVectorX* descriptor) {
  CHECK_NOTNULL(descriptor);
  CHECK_GT(pixel_samples_.size(), 0) << "You must call Initialize() first!";

  // Smooth the image using a gaussian kernal of 2.
  FloatImage blurred_image = image.AsGrayscaleImage();
  blurred_image.ApproximateGaussianBlur(2.0);

  return ComputeDescriptorInternal(blurred_image,
                                   keypoint,
                                   descriptor);
}

// Computes a descriptor at a single keypoint.
bool BriefDescriptorExtractor::ComputeDescriptors(
    const FloatImage& image,
    std::vector<Keypoint>* keypoints,
    std::vector<BinaryVectorX>* descriptors) {
  CHECK_NOTNULL(keypoints);
  CHECK_NOTNULL(descriptors)->clear();
  CHECK_GT(pixel_samples_.size(), 0) << "You must call Initialize() first!";

  // Smooth the image using a gaussian kernal of 2.
  FloatImage blurred_image = image.AsGrayscaleImage();
  blurred_image.ApproximateGaussianBlur(2.0);

  auto it = keypoints->begin();
  while (it != keypoints->end()) {
    BinaryVectorX descriptor;
    if (!ComputeDescriptorInternal(blurred_image,
                                   *it,
                                   &descriptor)) {
      it = keypoints->erase(it);
      continue;
    }

    descriptors->emplace_back(descriptor);
    ++it;
  }

  return descriptors->size() > 0;
}

bool BriefDescriptorExtractor::ComputeDescriptorInternal(
    const FloatImage& blurred_image,
    const Keypoint& keypoint,
    BinaryVectorX* descriptor) {
  // Do not consider keypoints that are too close to the border.
  if (KeypointIsTooCloseToBorder(keypoint,
                                 blurred_image.Rows(),
                                 blurred_image.Cols(),
                                 patch_sample_size_)) {
    return false;
  }

  // Extract the descriptor.
  const int size = (num_samples_ / 8) / sizeof(uint8_t);
  descriptor->resize(size);

  for (int i = 0; i < num_samples_; i += 128) {
    const int current_index = (i / 8) / sizeof(uint8_t);
    std::bitset<128>* bits =
        reinterpret_cast<std::bitset<128>*>(descriptor->data() + current_index);
    for (int j = 0; j < 128; j++) {
      const double pixel1 = blurred_image(
          keypoint.x() + pixel_samples_[i + j].pixel_points[0],
          keypoint.y() + pixel_samples_[i + j].pixel_points[1]);
      const double pixel2 = blurred_image(
          keypoint.x() + pixel_samples_[i + j].pixel_points[2],
          keypoint.y() + pixel_samples_[i + j].pixel_points[3]);
      (*bits)[j] = pixel1 < pixel2;
    }
  }
  return true;
}

bool BriefDescriptorExtractor::DetectAndExtractDescriptors(
    const FloatImage& image,
    std::vector<Keypoint>* keypoints,
    std::vector<BinaryVectorX>* descriptors) {
  BriskDetector detector;
  if (!detector.DetectKeypoints(image, keypoints)) {
    return false;
  }

  return this->ComputeDescriptors(image, keypoints, descriptors);
}

const std::vector<BriefSamplePair>&
BriefDescriptorExtractor::SamplePairs() const {
  return pixel_samples_;
}

}  // namespace theia
