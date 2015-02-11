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

#include <gflags/gflags.h>
#include <glog/logging.h>
#include <string>
#include "gtest/gtest.h"

#include "theia/image/descriptor/binary_descriptor.h"
#include "theia/image/image.h"
#include "theia/image/keypoint_detector/brisk_detector.h"
#include "theia/image/descriptor/freak_descriptor.h"
#include "theia/matching/distance.h"

DEFINE_string(test_img, "image/descriptor/img1.png",
              "Name of test image file.");
namespace theia {
namespace {
std::string img_filename = THEIA_DATA_DIR + std::string("/") + FLAGS_test_img;
}  // namespace

TEST(FreakDescriptor, Sanity) {
  FloatImage input_img(img_filename);

  // Get keypoints.
  BriskDetector brisk_detector;
  std::vector<Keypoint> brisk_keypoints;
  brisk_detector.DetectKeypoints(input_img, &brisk_keypoints);

  // For each keypoint, extract the freak descriptors.
  FreakDescriptorExtractor freak_extractor;
  EXPECT_TRUE(freak_extractor.Initialize());
  BinaryVectorX descriptor;
  // We need to make sure we pick a point away from the border so that a
  // descriptor can actually be extracted.
  CHECK_GT(brisk_keypoints.size(), 100);
  EXPECT_TRUE(freak_extractor.ComputeDescriptor(input_img,
                                                brisk_keypoints[100],
                                                &descriptor));

  std::vector<BinaryVectorX> freak_descriptors;
  EXPECT_TRUE(freak_extractor.ComputeDescriptors(input_img,
                                                 &brisk_keypoints,
                                                 &freak_descriptors));
}

}  // namespace theia
