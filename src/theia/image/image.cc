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

#include "theia/image/image.h"

#include <Eigen/Core>
#include <glog/logging.h>
#include <OpenImageIO/imagebuf.h>
#include <OpenImageIO/imagebufalgo.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "theia/util/util.h"

namespace theia {

// Read from file.
FloatImage::FloatImage(const std::string& filename) {
  Read(filename);
}

FloatImage::FloatImage(const FloatImage& image_to_copy) {
  CHECK(image_.copy(image_to_copy.image_));
}

FloatImage::FloatImage(const int width, const int height, const int channels) {
  OpenImageIO::ImageSpec image_spec(width, height, channels,
                                    OpenImageIO::TypeDesc::FLOAT);
  image_.reset(image_spec);
}

FloatImage::FloatImage(const int width, const int height, const int channels,
                       float* buffer)
    : image_(OpenImageIO::ImageSpec(width, height, channels,
                                    OpenImageIO::TypeDesc::FLOAT),
             reinterpret_cast<void*>(buffer)) {}

FloatImage::FloatImage(const OpenImageIO::ImageBuf& image) {
  image_.copy(image);
}

FloatImage& FloatImage::operator=(const FloatImage& image2) {
  image_.copy(image2.image_);
  return *this;
}

OpenImageIO::ImageBuf& FloatImage::GetOpenImageIOImageBuf() {
  return image_;
}

const OpenImageIO::ImageBuf& FloatImage::GetOpenImageIOImageBuf() const {
  return image_;
}

int FloatImage::Rows() const { return Height(); }

int FloatImage::Cols() const { return Width(); }

int FloatImage::Width() const { return image_.spec().width; }

int FloatImage::Height() const { return image_.spec().height; }

int FloatImage::Channels() const { return image_.nchannels(); }

void FloatImage::SetXY(const int x,
                       const int y,
                       const int c,
                       const float value) {
  DCHECK_GE(x, 0);
  DCHECK_LT(x, Width());
  DCHECK_GE(y, 0);
  DCHECK_LT(y, Height());
  DCHECK_GE(c, 0);
  DCHECK_LT(c, Channels());

  // Set the ROI to be the precise pixel location in the correct channel.
  OpenImageIO::ImageBuf::Iterator<float> it(image_, x, y, 0);
  it[c] = value;
}

void FloatImage::SetXY(const int x, const int y, const Eigen::Vector3f& rgb) {
  DCHECK_GE(x, 0);
  DCHECK_LT(x, Width());
  DCHECK_GE(y, 0);
  DCHECK_LT(y, Height());

  image_.setpixel(x, y, 0, rgb.data());
}

float FloatImage::GetXY(const int x, const int y, const int c) const {
  DCHECK_GE(x, 0);
  DCHECK_LT(x, Width());
  DCHECK_GE(y, 0);
  DCHECK_LT(y, Height());
  DCHECK_GE(c, 0);
  DCHECK_LT(c, Channels());
  return image_.getchannel(x, y, 0, c);
}

Eigen::Vector3f FloatImage::GetXY(const int x, const int y) const {
  DCHECK_GE(x, 0);
  DCHECK_LT(x, Width());
  DCHECK_GE(y, 0);
  DCHECK_LT(y, Height());
  Eigen::Vector3f rgb;
  image_.getpixel(x, y, 0, rgb.data());
  return rgb;
}

void FloatImage::SetRowCol(const int row,
                           const int col,
                           const int channel,
                           const float value) {
  SetXY(col, row, channel, value);
}

void FloatImage::SetRowCol(const int row,
                           const int col,
                           const Eigen::Vector3f& rgb) {
  SetXY(col, row, rgb);
}

float FloatImage::GetRowCol(const int row,
                            const int col,
                            const int channel) const {
  return GetXY(col, row, channel);
}

Eigen::Vector3f FloatImage::GetRowCol(const int row,
                                      const int col) const {
  return GetXY(col, row);
}

float FloatImage::BilinearInterpolate(const double x,
                                      const double y,
                                      const int c) const {
  DCHECK_GE(c, 0);
  DCHECK_LT(c, Channels());
  // The caller has to ensure that "pixel" has sufficient memory to store all
  // interpolated channels so we have to initialize it as an array.
  std::vector<float> pixel(Channels());
  image_.interppixel(x, y, pixel.data());
  return pixel[c];
}

Eigen::Vector3f FloatImage::BilinearInterpolate(const double x,
                                                const double y) const {
  Eigen::Vector3f pixel;
  image_.interppixel(x, y, pixel.data());
  return pixel;
}

void FloatImage::ConvertToGrayscaleImage() {
  if (Channels() == 1) {
    VLOG(2) << "Image is already a grayscale image. No conversion necessary.";
    return;
  }

  // Compute luminance via a weighted sum of R,G,B (assuming Rec709 primaries
  // and a linear scale)
  const float luma_weights[3] = {.2126, .7152, .0722};
  OpenImageIO::ImageBuf source = image_;
  image_.clear();
  OpenImageIO::ImageBufAlgo::channel_sum(image_, source, luma_weights);
}

void FloatImage::ConvertToRGBImage() {
  if (Channels() == 3) {
    VLOG(2) << "Image is already an RGB image. No conversion necessary.";
    return;
  }

  // Copy the single grayscale channel into r, g, and b.
  const OpenImageIO::ImageBuf source(image_);
  OpenImageIO::ImageSpec image_spec(Width(), Height(), 3,
                                    OpenImageIO::TypeDesc::FLOAT);
  image_.reset(image_spec);
  OpenImageIO::ImageBufAlgo::paste(image_, 0, 0, 0, 0, source);
  OpenImageIO::ImageBufAlgo::paste(image_, 0, 0, 0, 1, source);
  OpenImageIO::ImageBufAlgo::paste(image_, 0, 0, 0, 2, source);
}

FloatImage FloatImage::AsGrayscaleImage() const {
  if (Channels() == 1) {
    VLOG(2) << "Image is already a grayscale image. No conversion necessary.";
    return *this;
  }
  FloatImage gray_image(*this);
  gray_image.ConvertToGrayscaleImage();
  return gray_image;
}

FloatImage FloatImage::AsRGBImage() const {
  if (Channels() == 3) {
    VLOG(2) << "Image is already an RGB image. No conversion necessary.";
    return *this;
  }

  FloatImage rgb_image(*this);
  rgb_image.ConvertToRGBImage();
  return rgb_image;
}

void FloatImage::ScalePixels(float scale) {
  OpenImageIO::ImageBufAlgo::mul(image_, image_, scale);
}

void FloatImage::Read(const std::string& filename) {
  image_.reset(filename);
  image_.read(0, 0, true, OpenImageIO::TypeDesc::FLOAT);
}

void FloatImage::Write(const std::string& filename) const {
  image_.write(filename);
}

float* FloatImage::Data() {
  return reinterpret_cast<float*>(image_.localpixels());
}
const float* FloatImage::Data() const {
  return reinterpret_cast<const float*>(image_.localpixels());
}

FloatImage FloatImage::ComputeGradientX() const {
  float sobel_filter_x[9] = {-.125, 0, .125, -.25, 0, .25, -.125, 0, .125};
  OpenImageIO::ImageSpec spec(3, 3, 1, OpenImageIO::TypeDesc::FLOAT);
  OpenImageIO::ImageBuf kernel_x(spec, sobel_filter_x);
  OpenImageIO::ImageBuf gradient_x;
  OpenImageIO::ImageBufAlgo::convolve(gradient_x, image_, kernel_x, false);
  return FloatImage(gradient_x);
}

FloatImage FloatImage::ComputeGradientY() const {
  float sobel_filter_y[9] = {-.125, -.25, -.125, 0, 0, 0, .125, .25, .125};
  OpenImageIO::ImageSpec spec(3, 3, 1, OpenImageIO::TypeDesc::FLOAT);
  OpenImageIO::ImageBuf kernel_y(spec, sobel_filter_y);
  OpenImageIO::ImageBuf gradient_y;
  OpenImageIO::ImageBufAlgo::convolve(gradient_y, image_, kernel_y, false);
  return FloatImage(gradient_y);
}

FloatImage FloatImage::ComputeGradient() const {
  // Get Dx and Dy.
  float sobel_filter_x[9] = {-.125, 0, .125, -.25, 0, .25, -.125, 0, .125};
  float sobel_filter_y[9] = {-.125, -.25, -.125, 0, 0, 0, .125, .25, .125};

  OpenImageIO::ImageSpec spec(3, 3, 1, OpenImageIO::TypeDesc::FLOAT);
  OpenImageIO::ImageBuf kernel_x(spec, sobel_filter_x);
  OpenImageIO::ImageBuf kernel_y(spec, sobel_filter_y);

  OpenImageIO::ImageBuf gradient, gradient_x, gradient_y;
  OpenImageIO::ImageBufAlgo::convolve(gradient_x, image_, kernel_x, false);
  OpenImageIO::ImageBufAlgo::abs(gradient_x, gradient_x);
  OpenImageIO::ImageBufAlgo::convolve(gradient_y, image_, kernel_y, false);
  OpenImageIO::ImageBufAlgo::abs(gradient_y, gradient_y);
  OpenImageIO::ImageBufAlgo::add(gradient, gradient_x, gradient_y);

  return FloatImage(gradient);
}

void FloatImage::ApproximateGaussianBlur(const int kernel_size) {
  OpenImageIO::ImageBuf kernel;
  OpenImageIO::ImageBufAlgo::make_kernel(kernel, "gaussian",
                                         static_cast<float>(kernel_size),
                                         static_cast<float>(kernel_size));
  OpenImageIO::ImageBufAlgo::convolve(image_, image_, kernel);
}

void FloatImage::MedianFilter(const int patch_width) {
  CHECK(OpenImageIO::ImageBufAlgo::median_filter(image_, image_, patch_width));
}

void FloatImage::Integrate(FloatImage* integral) const {
  integral->ResizeRowsCols(Rows() + 1, Cols() + 1);
  for (int i = 0; i < Channels(); i++) {
    // Fill the first row with zeros.
    for (int x = 0; x < Width(); x++) {
      integral->SetXY(x, 0, i, 0);
    }
    for (int y = 1; y <= Height(); y++) {
      // This variable is to correct floating point round off.
      float sum = 0;
      integral->SetXY(0, y, i, 0);
      for (int x = 1; x <= Width(); x++) {
        sum += this->GetXY(x - 1, y - 1, i);
        integral->SetXY(x, y, i, integral->GetXY(x, y - 1, i) + sum);
      }
    }
  }
}

void FloatImage::Resize(int new_width, int new_height, int num_channels) {
  // If the image has not been initialized then initialize it with the image
  // spec. Otherwise resize the image and interpolate pixels accordingly.
  if (!image_.initialized()) {
    OpenImageIO::ImageSpec image_spec(new_width, new_height, num_channels,
                                      OpenImageIO::TypeDesc::FLOAT);
    image_.reset(image_spec);
  } else {
    OpenImageIO::ROI roi(0, new_width, 0, new_height, 0, 1, 0, num_channels);
    OpenImageIO::ImageBuf dst;
    CHECK(OpenImageIO::ImageBufAlgo::resize(dst, image_, nullptr, roi))
      << OpenImageIO::geterror();
    image_.copy(dst);
  }
}

void FloatImage::Resize(int new_width, int new_height) {
  Resize(new_width, new_height, Channels());
}

void FloatImage::ResizeRowsCols(int new_rows, int new_cols) {
  Resize(new_cols, new_rows, Channels());
}

void FloatImage::Resize(double scale) {
  Resize(static_cast<int>(scale * Width()),
         static_cast<int>(scale * Height()),
         Channels());
}

}  // namespace theia
