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

#include <gflags/gflags.h>
#include <glog/logging.h>
#include <bitset>
#include <string>
#include "gtest/gtest.h"

#include "theia/image/descriptor/binary_descriptor.h"
#include "theia/image/image.h"
#include "theia/image/keypoint_detector/brisk_detector.h"
#include "theia/image/descriptor/brief_descriptor.h"
#include "theia/matching/distance.h"

DEFINE_string(test_img, "image/descriptor/img1.png",
              "Name of test image file.");

namespace theia {
namespace {
std::string img_filename = THEIA_DATA_DIR + std::string("/") + FLAGS_test_img;
}  // namespace

TEST(BriefDescriptor, Sanity) {
  FloatImage input_img(img_filename);

  // Get keypoints.
  BriskDetector brisk_detector;
  std::vector<Keypoint> brief_keypoints;
  brisk_detector.DetectKeypoints(input_img, &brief_keypoints);

  // For each keypoint, extract the brief descriptors.
  BriefDescriptorExtractor brief_extractor(48, 16);
  brief_extractor.Initialize();
  std::vector<BinaryVectorX> brief_descriptors;
  EXPECT_TRUE(brief_extractor.ComputeDescriptors(input_img,
                                                 &brief_keypoints,
                                                 &brief_descriptors));

  EXPECT_GT(brief_descriptors.size(), 0);

  FloatImage blurred_image = input_img.AsGrayscaleImage();
  blurred_image.ApproximateGaussianBlur(2.0);
  // Check the value of each point.
  const auto& pixel_samples = brief_extractor.SamplePairs();
  for (int i = 0; i < brief_keypoints.size(); i++) {
    const auto& descriptor = brief_descriptors[i];

    const std::bitset<128>* brief_bitset =
        reinterpret_cast<const std::bitset<128>*>(descriptor.data());

    std::bitset<128> expected_bitset;
    for (int j = 0; j < pixel_samples.size(); j++) {
      const double pixel1 = blurred_image(
          brief_keypoints[i].x() + pixel_samples[j].pixel_points[0],
          brief_keypoints[i].y() + pixel_samples[j].pixel_points[1]);
      const double pixel2 = blurred_image(
          brief_keypoints[i].x() + pixel_samples[j].pixel_points[2],
          brief_keypoints[i].y() + pixel_samples[j].pixel_points[3]);
      expected_bitset[j] = pixel1 < pixel2;
    }

    EXPECT_EQ(*brief_bitset, expected_bitset);
  }
}

TEST(BriefDescriptor, TooCloseToBorder) {
  FloatImage input_img(img_filename);

  // Create keypoints that are too close to the border.
  std::vector<Keypoint> brief_keypoints(4);
  brief_keypoints[0].set_x(10);
  brief_keypoints[0].set_y(10);
  brief_keypoints[1].set_x(10);
  brief_keypoints[1].set_y(input_img.Height() - 10);
  brief_keypoints[2].set_x(input_img.Width() - 10);
  brief_keypoints[2].set_y(10);
  brief_keypoints[3].set_x(input_img.Width() - 10);
  brief_keypoints[3].set_y(input_img.Height() - 10);

  // For each keypoint, extract the brief descriptors.
  BriefDescriptorExtractor brief_extractor;
  brief_extractor.Initialize();
  std::vector<BinaryVectorX> brief_descriptors;
  EXPECT_FALSE(brief_extractor.ComputeDescriptors(input_img,
                                                  &brief_keypoints,
                                                  &brief_descriptors));
}

}  // namespace theia
