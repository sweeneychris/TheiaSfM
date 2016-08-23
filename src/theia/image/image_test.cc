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

#include <Eigen/Core>
#include <gflags/gflags.h>
#include <OpenImageIO/imagebuf.h>
#include <OpenImageIO/imagebufalgo.h>
#include <stdio.h>
#include <string>

#include "gtest/gtest.h"
#include "theia/image/image.h"
#include "theia/util/random.h"

DEFINE_string(test_img, "image/test1.jpg", "Name of test image file.");

namespace theia {
namespace {

RandomNumberGenerator rng(51);

std::string img_filename = THEIA_DATA_DIR + std::string("/") + FLAGS_test_img;

#define ASSERT_IMG_EQ(oiio_img, theia_img, rows, cols)      \
  oiio_img.read(0, 0, true, OpenImageIO::TypeDesc::FLOAT);  \
  ASSERT_EQ(oiio_img.oriented_width(), theia_img.Cols());   \
  ASSERT_EQ(oiio_img.oriented_height(), theia_img.Rows());  \
  ASSERT_EQ(oiio_img.nchannels(), theia_img.Channels());    \
  OpenImageIO::ImageBuf::ConstIterator<float> it(oiio_img); \
  for (; !it.done(); ++it) {                                \
    for (int c = 0; c < oiio_img.nchannels(); c++) {        \
      ASSERT_EQ(it[c], theia_img.GetXY(it.x(), it.y(), c)); \
    }                                                       \
  }

float Interpolate(const FloatImage& image,
                  const double x,
                  const double y,
                  const int c) {
  const float x_fix = x - 0.5;
  const float y_fix = y - 0.5;

  float intpart;
  const float s = std::modf(x_fix, &intpart);
  const int left = static_cast<int>(intpart);
  const float t = std::modf(y_fix, &intpart);
  const int top = static_cast<int>(intpart);

  const float v0 = image.GetXY(left, top, 0);
  const float v1 = image.GetXY(left + 1, top, 0);
  const float v2 = image.GetXY(left, top + 1, 0);
  const float v3 = image.GetXY(left + 1, top + 1, 0);

  return (1.0 - t) * (v0 * (1.0 - s) + v1 * s) + t * (v2 * (1.0 - s) + v3 * s);
}

}  // namespace

// Test that inputting the old fashioned way is the same as through our class.
TEST(Image, RGBInput) {
  OpenImageIO::ImageBuf oiio_img(img_filename.c_str());
  oiio_img.read();
  FloatImage theia_img(img_filename);

  int rows = oiio_img.oriented_height();
  int cols = oiio_img.oriented_width();

  // Assert each pixel value is exactly the same!
  ASSERT_IMG_EQ(oiio_img, theia_img, rows, cols);
}

// Test that width and height methods work.
TEST(Image, RGBColsRows) {
  OpenImageIO::ImageBuf oiio_img(img_filename.c_str());
  FloatImage theia_img(img_filename);

  int true_height = oiio_img.oriented_height();
  int true_width = oiio_img.oriented_width();

  ASSERT_EQ(theia_img.Cols(), true_width);
  ASSERT_EQ(theia_img.Rows(), true_height);
}

// Test that inputting the old fashioned way is the same as through our class.
TEST(Image, ConvertToGrayscaleImage) {
  OpenImageIO::ImageBuf oiio_img(img_filename.c_str());
  OpenImageIO::ImageBuf gray_img;
  const float luma_weights[3] = {.2126, .7152, .0722};
  OpenImageIO::ImageBufAlgo::channel_sum(gray_img, oiio_img, luma_weights);

  FloatImage theia_img(img_filename);
  theia_img.ConvertToGrayscaleImage();
  ASSERT_EQ(theia_img.Channels(), 1);

  int rows = oiio_img.oriented_height();
  int cols = oiio_img.oriented_width();

  // Assert each pixel value is exactly the same!
  ASSERT_IMG_EQ(gray_img, theia_img, rows, cols);
}

TEST(Image, ConvertToRGBImage) {
  OpenImageIO::ImageBuf oiio_img(img_filename.c_str());
  OpenImageIO::ImageBuf gray_img;
  const float luma_weights[3] = {.2126, .7152, .0722};
  OpenImageIO::ImageBufAlgo::channel_sum(gray_img, oiio_img, luma_weights);

  // This should result in an image with the grayscale image copied in each
  // channel.
  FloatImage rgb_img(img_filename.c_str());
  rgb_img.ConvertToGrayscaleImage();
  rgb_img.ConvertToRGBImage();

  ASSERT_EQ(rgb_img.Width(), gray_img.oriented_width());
  ASSERT_EQ(rgb_img.Height(), gray_img.oriented_height());
  ASSERT_EQ(rgb_img.Channels(), 3);

  // Check that all channels have equal value and that the value is equal to the
  // grayscale image.
  for (OpenImageIO::ImageBuf::ConstIterator<float> it(gray_img);
       !it.done();
       ++it) {
    ASSERT_EQ(it[0], rgb_img.GetXY(it.x(), it.y(), 0));
    ASSERT_EQ(it[0], rgb_img.GetXY(it.x(), it.y(), 1));
    ASSERT_EQ(it[0], rgb_img.GetXY(it.x(), it.y(), 2));
  }
}

TEST(Image, BillinearInterpolate) {
  static const int kNumTrials = 10;
  static const float kTolerance = 1e-2;

  FloatImage theia_img(img_filename);
  theia_img.ConvertToGrayscaleImage();
  for (int i = 0; i < kNumTrials; i++) {
    const double x = rng.RandDouble(1.0, theia_img.Width() - 2);
    const double y = rng.RandDouble(1.0, theia_img.Height() - 2);
    const float pixel = Interpolate(theia_img, x, y, 0);
    const float pixel2 = theia_img.BilinearInterpolate(x, y, 0);
    EXPECT_NEAR(pixel, pixel2, kTolerance);
  }
}

TEST(Image, ScalePixels) {
  static const float kTolerance = 1e-2;
  static const float kScaleFactor = 1.1;

  FloatImage theia_img(img_filename);
  FloatImage scaled_img(img_filename);
  scaled_img.ScalePixels(kScaleFactor);

  for (int y = 0; y < theia_img.Height(); y++) {
    for (int x = 0; x < theia_img.Width(); x++) {
      ASSERT_EQ(kScaleFactor * theia_img.GetXY(x, y, 0),
                scaled_img.GetXY(x, y, 0));
    }
  }
}

TEST(Image, Resize) {
  static const int kWidth = 800;
  static const int kHeight = 600;

  FloatImage theia_img(img_filename);
  theia_img.Resize(kWidth, kHeight);

  // Make sure the image was resized appropriately.
  EXPECT_EQ(theia_img.Width(), kWidth);
  EXPECT_EQ(theia_img.Height(), kHeight);
}

TEST(Image, ResizeUninitialized) {
  static const int kWidth = 800;
  static const int kHeight = 600;

  FloatImage theia_img;
  theia_img.Resize(kWidth, kHeight);

  // Make sure the image was resized appropriately.
  EXPECT_EQ(theia_img.Width(), kWidth);
  EXPECT_EQ(theia_img.Height(), kHeight);

  // Make sure the resizing still works when converting an uninitialized image
  // to RGB first.
  FloatImage theia_img2;
  theia_img2.ConvertToRGBImage();
  theia_img2.Resize(kWidth, kHeight);

  // Make sure the image was resized appropriately.
  EXPECT_EQ(theia_img2.Width(), kWidth);
  EXPECT_EQ(theia_img2.Height(), kHeight);
  EXPECT_EQ(theia_img2.Channels(), 3);
}

}  // namespace theia
