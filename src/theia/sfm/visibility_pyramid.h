// Copyright (C) 2017 The Regents of the University of California (Regents).
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
// Author: Chris Sweeney (sweeney.chris.m@gmail.com)

#ifndef THEIA_SFM_VISIBILITY_PYRMAID_H_
#define THEIA_SFM_VISIBILITY_PYRAMID_H_

#include <Eigen/Core>
#include <vector>

namespace theia {

// An implementation of the visibility pyramid used for "next best view" ranking
// as described in: "Structure-from-Motion Revisited" by Johannes Schonberger
// and Jan-Michael Frahm (CVPR 2016). The main idea is that well-constrained
// views have a large number of observations and good spatial distribution of
// those observations. The VisibilityPyramid is used to compute a score for a
// view based on the aforementioned criteria. The view is divided into grid
// cells and a score is assigned to the grid based on the number of occupied
// cells. The score is accumulated over a pyrammid to capture coarse and
// fine-granded distribution statistics.
//
// An online version of this class is available with the COLMAP library, and
// this implementation was slightly inspired by that:
//   https://github.com/colmap/colmap/blob/master/src/base/visibility_pyramid.h
class VisibilityPyramid {
 public:
  // The inputs are the view/image width and height, as well as the number of
  // desired levels in the image pyramid.
  VisibilityPyramid(const int width,
                    const int height,
                    const int num_pyramid_levels);

  // Add a point to the visibility pyramid.
  void AddPoint(const Eigen::Vector2d& point);

  // Compute the score of the visibility pyramid. Higher scores indicate that
  // the view is better constrained by the points.
  int ComputeScore() const;

 private:
  const int width_, height_, num_pyramid_levels_, max_cells_in_dimension_;
  // The pyramid represents all levels of image grids that keep track of the
  // features. The score of the pyramid may be efficiently computed with the
  // count() method that Eigen provides to determine occupancy.
  //
  // The pyramid is stored from coarse to fine and is indexed as (x, y).
  std::vector<Eigen::MatrixXi> pyramid_;
};
}  // namespace theia

#endif  // THEIA_SFM_VISIBILITY_PYRAMID_H_
