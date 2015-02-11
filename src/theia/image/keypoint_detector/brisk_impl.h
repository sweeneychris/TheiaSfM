/*
  Brisk reference implementation, modified to remove OpenCV dependency and use
  Theia's Image class instead.
  Changes Copyright (C) 2013 Chris Sweeney (cmsweeney@cs.ucsb.edu)

  BRISK - Binary Robust Invariant Scalable Keypoints
  Reference implementation of
  [1] Stefan Leutenegger,Margarita Chli and Roland Siegwart, BRISK:
  Binary Robust Invariant Scalable Keypoints, in Proceedings of
  the IEEE International Conference on Computer Vision (ICCV2011).

  Copyright (C) 2011  The Autonomous Systems Lab (ASL), ETH Zurich,
  Stefan Leutenegger, Simon Lynen and Margarita Chli.

  This file is part of BRISK.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are met:
  * Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.
  * Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.
  * Neither the name of the ASL nor the names of its contributors may be
  used to endorse or promote products derived from this software without
  specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
  AND
  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS BE LIABLE FOR ANY
  DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef THEIA_IMAGE_KEYPOINT_DETECTOR_BRISK_IMPL_H_
#define THEIA_IMAGE_KEYPOINT_DETECTOR_BRISK_IMPL_H_

#include <agast/oast9_16.h>
#include <agast/agast7_12s.h>
#include <agast/agast5_8.h>
#include <agast/cvWrapper.h>
#include <Eigen/Core>
#include <vector>

#include "theia/image/image.h"
#include "theia/image/keypoint_detector/keypoint.h"

#ifndef M_PI
#define M_PI 3.141592653589793
#endif

namespace theia {

// a layer in the Brisk detector pyramid
class BriskLayer {
 public:
  // constructor arguments
  struct CommonParams {
    static const int HALFSAMPLE = 0;
    static const int TWOTHIRDSAMPLE = 1;
  };
  // construct a base layer
  BriskLayer(const Image<unsigned char>& img, float scale = 1.0f,
             float offset = 0.0f);
  // derive a layer
  BriskLayer(const BriskLayer& layer, int mode);

  // Fast/Agast without non-max suppression
  void getAgastPoints(uint8_t threshold,
                      std::vector<OpenCVPoint>& keypoints); // NOLINT

  // get scores - NOTE: this is in layer coordinates, not scale=1 coordinates!
  inline uint8_t getAgastScore(int x, int y, uint8_t threshold);
  inline uint8_t getAgastScore(float xf, float yf, uint8_t threshold,
                               float scale = 1.0f);
  inline uint8_t getAgastScore_5_8(int x, int y, uint8_t threshold);

  // accessors
  inline const Image<unsigned char>& img() const { return img_; }
  inline const Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic>&
  scores() const {
    return scores_;
  }

  inline float scale() const { return scale_; }
  inline float offset() const { return offset_; }

 private:
  // access gray values (smoothed/interpolated)
  inline uint8_t value(
      const Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic>& mat,
      float xf, float yf, float scale);
  // the image
  Image<unsigned char> img_;
  // its Fast scores
  Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> scores_;
  // coordinate transformation
  float scale_;
  float offset_;
  // agast
  std::shared_ptr<agast::OastDetector9_16> oastDetector_;
  std::shared_ptr<agast::AgastDetector5_8> agastDetector_5_8_;
};

class BriskScaleSpace {
 public:
  // construct telling the octaves number:
  explicit BriskScaleSpace(uint8_t _octaves = 3);
  ~BriskScaleSpace() {}

  // construct the image pyramids
  void constructPyramid(const Image<unsigned char>& image);

  // get Keypoints
  void getKeypoints(const uint8_t _threshold,
                    std::vector<Keypoint>* keypoints);

 private:
  // nonmax suppression:
  inline bool isMax2D(const uint8_t layer, const int x_layer,
                      const int y_layer);
  // 1D (scale axis) refinement:
  // around octave
  inline float refine1D(const float s_05, const float s0, const float s05,
                        float& max);  // NOLINT
  // around intra
  inline float refine1D_1(const float s_05, const float s0, const float s05,
                          float& max);  // NOLINT
  // around octave 0 only
  inline float refine1D_2(const float s_05, const float s0, const float s05,
                          float& max);  // NOLINT
  // 2D maximum refinement:
  inline float subpixel2D(const int s_0_0, const int s_0_1, const int s_0_2,
                          const int s_1_0, const int s_1_1, const int s_1_2,
                          const int s_2_0, const int s_2_1, const int s_2_2,
                          float& delta_x, float& delta_y);  // NOLINT
  // 3D maximum refinement centered around (x_layer,y_layer)
  inline float refine3D(const uint8_t layer, const int x_layer,
                        const int y_layer,
                        float& x, float& y,  // NOLINT
                        float& scale, bool& ismax);  // NOLINT

  // interpolated score access with recalculation when needed:
  inline int getScoreAbove(const uint8_t layer, const int x_layer,
                           const int y_layer);
  inline int getScoreBelow(const uint8_t layer, const int x_layer,
                           const int y_layer);

  // return the maximum of score patches above or below
  inline float getScoreMaxAbove(const uint8_t layer, const int x_layer,
                                const int y_layer, const int threshold,
                                bool& ismax, float& dx, float& dy);  // NOLINT
  inline float getScoreMaxBelow(const uint8_t layer, const int x_layer,
                                const int y_layer, const int threshold,
                                bool& ismax, float& dx, float& dy);  // NOLINT

  // the image pyramids:
  uint8_t layers_;
  std::vector<BriskLayer> pyramid_;

  // Agast:
  uint8_t threshold_;
  uint8_t safeThreshold_;

  // some constant parameters:
  static const float safetyFactor_;
  static const float basicSize_;
};
}  // namespace theia

#endif  // THEIA_IMAGE_KEYPOINT_DETECTOR_BRISK_IMPL_H_
