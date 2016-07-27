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

#include "gtest/gtest.h"

#include "theia/sfm/camera_intrinsics_prior.h"
#include "theia/sfm/exif_reader.h"

namespace theia {

std::string exif_img_filename = THEIA_DATA_DIR + std::string("/image/exif.jpg");

TEST(ExtractEXIFMetadata, FocalLengthEXIF) {
  CameraIntrinsicsPrior camera_intrinsics_prior;
  ExifReader exif_reader;
  EXPECT_TRUE(exif_reader.ExtractEXIFMetadata(exif_img_filename,
                                              &camera_intrinsics_prior));
  EXPECT_TRUE(camera_intrinsics_prior.focal_length.is_set);
  EXPECT_NEAR(camera_intrinsics_prior.focal_length.value[0], 1304.84, 0.1);
  EXPECT_EQ(camera_intrinsics_prior.image_width, 960);
  EXPECT_EQ(camera_intrinsics_prior.image_height, 1280);
}

}  // namespace theia
