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

#ifndef THEIA_IMAGE_IMAGE_H_
#define THEIA_IMAGE_IMAGE_H_

#include <Eigen/Core>
#include <glog/logging.h>
#include <OpenImageIO/imagebuf.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "theia/util/util.h"

namespace theia {

// A basic wrapper class for handling images. The images are always converted to
// floating point type with pixel values ranging from 0 to 1.0. The number of
// channels is dynamically controlled and methods that access pixels requires
// the caller to choose to appropriate RGB vs grayscale method.
class FloatImage {
 public:
  FloatImage() {}

  // Read from file.
  explicit FloatImage(const std::string& filename);
  FloatImage(const int width, const int height, const int channels);

  // Copy function. This is a deep copy of the image.
  FloatImage(const FloatImage& image_to_copy);

  explicit FloatImage(const OpenImageIO::ImageBuf& image);

  ~FloatImage() {}

  // Image information
  int Rows() const;
  int Cols() const;
  int Width() const;
  int Height() const;
  int Channels() const;

  // Set the pixel color value at (x, y) in channel c.
  void SetXY(const int x,
             const int y,
             const int c,
             const float value);
  // Set the rgb color values at the pixel (x, y). This assumes (with no
  // checks!) that the image is an rgb image.
  void SetXY(const int x,
             const int y,
             const Eigen::Vector3f& rgb);

  // Get the pixel value at the given (x, y) position and channel.
  float GetXY(const int x, const int y, const int channel) const;
  // Get the RGB value at the given pixel. This assumes (with no checks!) that
  // the image is indeed an RGB image.
  Eigen::Vector3f GetXY(const int x, const int y) const;

  // Set the pixel value at row and column in channel c.
  void SetRowCol(const int row,
                 const int col,
                 const int channel,
                 const float value);
  // Set the rgb color values at the given row and column. This assumes (with no
  // checks!) that the image is an rgb image.
  void SetRowCol(const int row,
                 const int col,
                 const Eigen::Vector3f& rgb);
  // Get the pixel value at the given location and channel.
  float GetRowCol(const int row, const int col, const int channel) const;

  // Get the RGB value at the given pixel. This assumes (with no checks!) that
  // the image is indeed an RGB image.
  Eigen::Vector3f GetRowCol(const int row, const int col) const;

  // Get the pixel value at a non-discrete location with bilinear interpolation.
  float BilinearInterpolate(const double x, const double y, const int c) const;

  // Get the pixel value at a non-discrete location with bilinear interpolation.
  Eigen::Vector3f BilinearInterpolate(const double x, const double y) const;

  // Convert to other image types.
  FloatImage AsGrayscaleImage() const;
  FloatImage AsRGBImage() const;
  void ConvertToGrayscaleImage();
  void ConvertToRGBImage();

  // Scale all the pixel values by a scale factor.
  void ScalePixels(float scale);

  // Write image to file.
  void Read(const std::string& filename);
  void Write(const std::string& filename) const;

  // Get a pointer to the data.
  float* Data();
  const float* Data() const;

  // Computes the gradient in x and y and returns the summation to obtain the
  // gradient magnitude at each pixel.
  FloatImage ComputeGradient() const;

  // Compute the integral image where pixel (x, y) is equal to the sum of all
  // values in the rectangle from (0, 0) to (x, y) non-inclusive. This means
  // that the first row and column are all zeros, and the returned integral
  // image is one pixel wider and taller than the caller.
  //
  // NOTE: This method should be called with precise number types such as double
  // otherwise floating roundoff errors are sure to occur.
  void Integrate(FloatImage* integral) const;

  // Computes a fast approximate gaussian blur of te image.
  void ApproximateGaussianBlur(const double sigma);

  // Resize using a Lanczos 3 filter.
  void Resize(int new_width, int new_height);
  void ResizeRowsCols(int new_rows, int new_cols);
  void Resize(double scale);

 protected:
  //template<class AnyType> friend class Image;
  // friend class ImageCanvas;

  OpenImageIO::ImageBuf image_;
};
}  // namespace theia

#endif  // THEIA_IMAGE_IMAGE_H_
