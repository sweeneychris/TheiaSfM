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

#include <cimg/CImg.h>
#include <Eigen/Core>
#include <gflags/gflags.h>
#include <stdio.h>
#include <string>

#include "gtest/gtest.h"
#include "theia/image/image.h"
#include "theia/util/random.h"

DEFINE_string(test_img, "image/test1.jpg", "Name of test image file.");

namespace theia {
namespace {
using cimg_library::CImg;

std::string img_filename = THEIA_DATA_DIR + std::string("/") + FLAGS_test_img;

#define ASSERT_RGB_IMG_EQ(cimg_img, theia_img, rows, cols)     \
  ASSERT_EQ(cimg_img.width(), theia_img.Cols());               \
  ASSERT_EQ(cimg_img.height(), theia_img.Rows());              \
  ASSERT_EQ(cimg_img.spectrum(), theia_img.Channels());        \
  ASSERT_EQ(theia_img.Channels(), 3);                          \
  ASSERT_EQ(cimg_img.depth(), 1);                              \
  for (int i = 0; i < cols; i++) {                             \
    for (int j = 0; j < rows; j++) {                           \
      ASSERT_EQ(cimg_img(i, j, 0, 0), theia_img(i, j, 0));     \
      ASSERT_EQ(cimg_img(i, j, 0, 1), theia_img(i, j, 1));     \
      ASSERT_EQ(cimg_img(i, j, 0, 2), theia_img(i, j, 2));     \
    }                                                          \
  }

#define ASSERT_GRAY_IMG_EQ(cimg_img, theia_img, rows, cols)     \
  ASSERT_EQ(cimg_img.width(), theia_img.Cols());                \
  ASSERT_EQ(cimg_img.height(), theia_img.Rows());               \
  ASSERT_EQ(cimg_img.spectrum(), theia_img.Channels());         \
  ASSERT_EQ(theia_img.Channels(), 1);                           \
  ASSERT_EQ(cimg_img.depth(), 1);                               \
  for (int i = 0; i < cols; i++)                                \
    for (int j = 0; j < rows; j++)                              \
      ASSERT_EQ(cimg_img(i, j), theia_img(i, j));               \

float Interpolate(const FloatImage& image,
                  const double x,
                  const double y,
                  const int c) {
  const int left = std::floor(x);
  const int right = std::ceil(x);
  const int top = std::floor(y);
  const int bottom = std::ceil(y);
  return image(left, top, c) * (right - x) * (bottom - y) +
         image(left, bottom, c) * (right - x) * (y - top) +
         image(right, top, c) * (x - left) * (bottom - y) +
         image(right, bottom, c) * (x - left) * (y - top);
}

}  // namespace

// Test that inputting the old fashioned way is the same as through our class.
TEST(Image, RGBInput) {
  CImg<float> cimg_img(img_filename.c_str());
  FloatImage theia_img(img_filename);

  int rows = cimg_img.height();
  int cols = cimg_img.width();

  // Assert each pixel value is exactly the same!
  ASSERT_RGB_IMG_EQ(cimg_img, theia_img, rows, cols);
}

// Test that width and height methods work.
TEST(Image, RGBColsRows) {
  CImg<float> cimg_img(img_filename.c_str());
  FloatImage theia_img(img_filename);

  int true_height = cimg_img.height();
  int true_width = cimg_img.width();

  ASSERT_EQ(theia_img.Cols(), true_width);
  ASSERT_EQ(theia_img.Rows(), true_height);
}

// Test that inputting the old fashioned way is the same as through our class.
TEST(Image, ConvertToGrayscaleImage) {
  CImg<float> cimg_img(img_filename.c_str());
  CImg<float> gray_img(cimg_img.RGBtoYCbCr().channel(0));
  FloatImage theia_img(img_filename);
  theia_img.ConvertToGrayscaleImage();

  int rows = cimg_img.height();
  int cols = cimg_img.width();

  // Assert each pixel value is exactly the same!
  ASSERT_GRAY_IMG_EQ(gray_img, theia_img, rows, cols);
}

TEST(Image, ConvertToRGBImage) {
  const CImg<float> cimg_img(img_filename.c_str());
  const CImg<float>& gray_cimg = cimg_img.get_RGBtoYCbCr().get_channel(0);

  CImg<float> rgb_img = gray_cimg;
  rgb_img.resize(cimg_img.width(), cimg_img.height(), cimg_img.depth(),
                 3);

  cimg_forXY(rgb_img, x, y) {
    CHECK_EQ(rgb_img(x, y, 0, 0), rgb_img(x, y, 0, 1));
    CHECK_EQ(rgb_img(x, y, 0, 0), rgb_img(x, y, 0, 2));
  }

  FloatImage theia_img(img_filename);
  theia_img.ConvertToGrayscaleImage();
  theia_img.ConvertToRGBImage();

  int rows = cimg_img.height();
  int cols = cimg_img.width();

  // Assert each pixel value is exactly the same!
  ASSERT_RGB_IMG_EQ(rgb_img, theia_img, rows, cols);
}

TEST(Image, IntegralImage) {
  const FloatImage img = FloatImage(img_filename).AsGrayscaleImage();
  Image<double> integral_img;
  img.Integrate(&integral_img);

  // Check the integral image over 100 trials;
  for (int i = 0; i < 1000; i++) {
    InitRandomGenerator();
    const int x = RandInt(1, img.Cols());
    const int y = RandInt(1, img.Rows());

    // Check the integral.
    double sum = 0;
    for (int r = 0; r < y; r++) {
      for (int c = 0; c < x; c++) {
        sum += img(c, r);
      }
    }

    EXPECT_DOUBLE_EQ(integral_img(x, y), sum);
  }
}

TEST(Image, BillinearInterpolate) {
  static const int kNumTrials = 10;
  static const float kTolerance = 1e-2;

  FloatImage theia_img(img_filename);
  InitRandomGenerator();
  for (int i = 0; i < kNumTrials; i++) {
    const double x = RandDouble(1.0, theia_img.Width() - 2);
    const double y = RandDouble(1.0, theia_img.Height() - 2);
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
      EXPECT_DOUBLE_EQ(kScaleFactor * theia_img(x, y), scaled_img(x, y));
    }
  }
}


}  // namespace theia
