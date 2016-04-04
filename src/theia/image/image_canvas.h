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

#ifndef THEIA_IMAGE_IMAGE_CANVAS_H_
#define THEIA_IMAGE_IMAGE_CANVAS_H_

#include <Eigen/Core>

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "theia/image/image.h"
#include "theia/image/keypoint_detector/keypoint.h"
#include "theia/matching/feature_correspondence.h"
#include "theia/util/random.h"
#include "theia/util/util.h"

namespace theia {

// This class allows for drawing on top of one or many images for visualization
// without modifying the underlying image. This is useful for feature detection,
// descriptor extraction, matching, etc. Add images to the canvas using
// AddImage, then draw on the canvas appropriately. There are two versions of
// all functions: one where you can specify drawing coordinates relative to the
// image at image_index, and another where you can specify the raw canvas
// coordinates. The likely use case of the latter is when you have an image
// canvas with only one image.
//
// NOTE: ImageCanvas objects are always an rgb image underneath so that pixels
// drawn onto the canvas can be color, even if the image is not. Again, this can
// be useful for visualization and matching.
struct RGBPixel{
  float r;
  float g;
  float b;
  RGBPixel(const float r, const float g, const float b) : r(r), g(g), b(b) {}
  RGBPixel() {}
};

class ImageCanvas {
 public:
  ImageCanvas() {}
  ~ImageCanvas() {}

  // Add an image to the canvas such that all the images that have been added
  // are now side-by-side on the canvas. This is useful for feature matching.
  int AddImage(const FloatImage& image);

  // Draw a circle in the image at image_index.
  void DrawCircle(int image_index, int x, int y, int radius,
                  const RGBPixel& color);
  // Draw a circle onto the canvas.
  void DrawCircle(int x, int y, int radius, const RGBPixel& color);

  // Draw a line in the image at image_index.
  void DrawLine(int image_index, int x1, int y1, int x2, int y2,
                const RGBPixel& color);
  void DrawLine(int image_index1, int x1, int y1, int image_index2, int x2,
                int y2, const RGBPixel& color);
  // Draw a line onto the image canvas.
  void DrawLine(int x1, int y1, int x2, int y2, const RGBPixel& color);

  // Draw feature in the image at image_index. We template these methods so that
  // you can use Keypoint or Descriptor types.
  // Draw feature in the image at image_index.
  template <class Feature>
  void DrawFeature(int image_index, const Feature& feature,
                   const RGBPixel& color, double scale);

  // Draw the feature onto the image canvas.
  template<class Feature>
  void DrawFeature(const Feature& feature, const RGBPixel& color,
                   double scale = 10.0);

  // Draw feature in the image at image_index.
  template<class Feature>
  void DrawFeatures(int image_index, const std::vector<Feature>& features,
                    const std::vector<RGBPixel>& colors, double scale = 10.0);
  template<class Feature>
  void DrawFeatures(int image_index, const std::vector<Feature>& features,
                    const RGBPixel& color, double scale = 10.0);
  // Draw the feature onto the image canvas.
  template<class Feature>
  void DrawFeatures(const std::vector<Feature>& features,
                    const std::vector<RGBPixel>& colors, double scale = 10.0);
  template<class Feature>
  void DrawFeatures(const std::vector<Feature>& features,
                    const RGBPixel& color, double scale = 10.0);

  // Draw matching features in the image by drawing a line from features1[i]
  // to features2[i].
  void DrawMatchedFeatures(int image_index1,
                           int image_index2,
                           const std::vector<FeatureCorrespondence>& matches,
                           double scale = 10.0);

  // Write the image canvas to a file.
  void Write(const std::string& output_name);

 private:
  void DrawFeature(int image_index,
                   int x, int y,
                   int radius,
                   double orientation,
                   const RGBPixel& color);
  // The local copy of the canvas. This can be comprised of image(s), shape(s),
  // and more.
  cimg_library::CImg<float> image_;

  // Contains the starting x coordinate of the image corresponding to the
  // index. This makes it easy to draw points relative to a particular image
  // when there are multiple images on the canvas.
  std::vector<int> pixel_offsets_;

  DISALLOW_COPY_AND_ASSIGN(ImageCanvas);
};

// ------------------ Implementation of template funcs ------------------ //

template <>
void ImageCanvas::DrawFeature<Keypoint>(int image_index,
                                        const Keypoint& feature,
                                        const RGBPixel& color, double scale);

template <>
void ImageCanvas::DrawFeature<Eigen::Vector2d>(int image_index,
                                               const Eigen::Vector2d& feature,
                                               const RGBPixel& color,
                                               double scale);

template <class Feature>
void ImageCanvas::DrawFeature(int image_index, const Feature& feature,
                              const RGBPixel& color, double scale) {
  double radius = feature.has_strength() ? feature.strength() * scale : scale;
  double angle = feature.has_orientation() ? feature.orientation() : 0.0;
  DrawFeature(image_index, feature.x(), feature.y(), radius, angle, color);
}

// Draw the feature onto the image canvas.
template<class Feature>
void ImageCanvas::DrawFeature(const Feature& feature, const RGBPixel& color,
                              double scale) {
  DrawFeature(0, feature, color, scale);
}

// Draw feature in the image at image_index.
template<class Feature>
void ImageCanvas::DrawFeatures(int image_index,
                               const std::vector<Feature>& features,
                               const std::vector<RGBPixel>& colors,
                               double scale) {
  for (int i = 0; i < features.size(); i++) {
    if (features[i] != nullptr)
      DrawFeature(image_index, features[i], colors[i], scale);
  }
}

// Draw feature in the image at image_index.
template<class Feature>
void ImageCanvas::DrawFeatures(int image_index,
                               const std::vector<Feature>& features,
                               const RGBPixel& color,
                               double scale) {
  for (int i = 0; i < features.size(); i++) {
    DrawFeature(image_index, features[i], color, scale);
  }
}

// Draw the feature onto the image canvas.
template<class Feature>
void ImageCanvas::DrawFeatures(const std::vector<Feature>& features,
                               const std::vector<RGBPixel>& colors,
                               double scale) {
  DrawFeatures(0, features, colors, scale);
}

// Draw the feature onto the image canvas.
template<class Feature>
void ImageCanvas::DrawFeatures(const std::vector<Feature>& features,
                               const RGBPixel& color,
                               double scale) {
  DrawFeatures(0, features, color, scale);
}

}  // namespace theia

#endif  // THEIA_IMAGE_IMAGE_CANVAS_H_
