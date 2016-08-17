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

#include "theia/image/image_canvas.h"

#include <cimg/CImg.h>
#include <glog/logging.h>

#include <cmath>
#include <string>
#include <vector>

#include "theia/image/image.h"
#include "theia/image/keypoint_detector/keypoint.h"

#define _USE_MATH_DEFINES

namespace theia {

// Add an image to the canvas such that all the images that have been added
// are now side-by-side on the canvas. This is useful for feature matching.
int ImageCanvas::AddImage(const FloatImage& image) {
  FloatImage rgb_image = image.AsRGBImage();

  if (pixel_offsets_.size() == 0) {
    image_ = rgb_image.image_;
    pixel_offsets_.push_back(0);
  } else {
    pixel_offsets_.push_back(image_.width());
    image_.append(rgb_image.image_);
  }
  return pixel_offsets_.size() - 1;
}

// Draw a circle in the image at image_index.
void ImageCanvas::DrawCircle(int image_index, int x, int y, int radius,
                             const RGBPixel& color) {
  CHECK_GT(pixel_offsets_.size(), image_index)
      << "Trying to draw in an image index that does not exist!";
  DrawCircle(pixel_offsets_[image_index] + x, y, radius, color);
}
// Draw a circle onto the canvas.
void ImageCanvas::DrawCircle(int x, int y, int radius, const RGBPixel& color) {
  float cimg_color[3] = { color.r, color.g, color.b };
  image_.draw_circle(x, y, radius, cimg_color, 1.0, 0L);
}

// Draw a line in the image at image_index.
void ImageCanvas::DrawLine(int image_index, int x1, int y1, int x2, int y2,
                           const RGBPixel& color) {
  CHECK_GT(pixel_offsets_.size(), image_index)
      << "Trying to draw in an image index that does not exist!";
  DrawLine(pixel_offsets_[image_index] + x1, y1,
           pixel_offsets_[image_index] + x2, y2, color);
}

void ImageCanvas::DrawLine(int image_index1, int x1, int y1, int image_index2,
                           int x2, int y2, const RGBPixel& color) {
  DrawLine(pixel_offsets_[image_index1] + x1, y1,
           pixel_offsets_[image_index2] + x2, y2, color);
}

// Draw a line onto the image canvas.
void ImageCanvas::DrawLine(int x1, int y1, int x2, int y2,
                           const RGBPixel& color) {
  double cimg_color[3] = { color.r, color.g, color.b };
  image_.draw_line(x1, y1, x2, y2, cimg_color);
}

void ImageCanvas::DrawFeature(int image_index,
                              int x, int y,
                              int radius,
                              double orientation,
                              const RGBPixel& color) {
  CHECK_GT(pixel_offsets_.size(), image_index)
      << "Trying to draw in an image index that does not exist!";
  // Draw circle at keypoint with size scale*strength.
  DrawCircle(image_index, x, y, radius, color);

  // Draw line in direction of the orientation if applicable.
  DrawLine(image_index,
           x, y,
           x + radius * cos(orientation), y + radius * sin(orientation),
           color);
}

// Write the image canvas to a file.
void ImageCanvas::Write(const std::string& output_name) {
  image_.save(output_name.c_str());
}

template <>
void ImageCanvas::DrawFeature(int image_index, const Keypoint& feature,
                              const RGBPixel& color, double scale) {
  double radius = feature.has_strength() ? feature.strength() * scale : scale;
  double angle = feature.has_orientation() ? feature.orientation() : 0.0;
  DrawFeature(image_index, feature.x(), feature.y(), radius, angle, color);
}

template <>
void ImageCanvas::DrawFeature<Eigen::Vector2d>(int image_index,
                                               const Eigen::Vector2d& feature,
                                               const RGBPixel& color,
                                               double scale) {
  double radius = scale;
  double angle = 0.0;
  DrawFeature(image_index, feature.x(), feature.y(), radius, angle, color);
}

// Draw matching features in the image.
void ImageCanvas::DrawMatchedFeatures(
    int image_index1,
    int image_index2,
    const std::vector<FeatureCorrespondence>& matches,
    double scale) {
  CHECK_GT(pixel_offsets_.size(), std::max(image_index1, image_index2));

  RandomNumberGenerator rng;
  for (int i = 0; i < matches.size(); i++) {
    const Feature& base = matches[i].feature1;
    const Feature& match = matches[i].feature2;
    RGBPixel color(rng.RandDouble(0, 255.0),
                   rng.RandDouble(0, 255.0),
                   rng.RandDouble(0, 255.0));
    DrawFeature(image_index1, base, color, scale);
    DrawFeature(image_index2, match, color, scale);
    DrawLine(pixel_offsets_[image_index1] + base.x(),
             base.y(),
             pixel_offsets_[image_index2] + match.x(),
             match.y(),
             color);
  }
}

}  // namespace theia
