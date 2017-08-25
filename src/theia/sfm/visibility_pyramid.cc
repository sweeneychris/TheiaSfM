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

#include "theia/sfm/visibility_pyramid.h"

#include <Eigen/Core>
#include <glog/logging.h>
#include <vector>

#include "theia/math/util.h"

namespace theia {

// The inputs are the view/image width and height, as well as the number of
// desired levels in the image pyramid.
VisibilityPyramid::VisibilityPyramid(const int width,
                                     const int height,
                                     const int num_pyramid_levels)
    : width_(width),
      height_(height),
      num_pyramid_levels_(num_pyramid_levels),
      max_cells_in_dimension_(1 << num_pyramid_levels) {
  CHECK_GT(width_, 0);
  CHECK_GT(height_, 0);
  CHECK_GT(num_pyramid_levels_, 0);

  pyramid_.resize(num_pyramid_levels);
  // Create the occupancy pyramid such that the coarsest level is a 2x2 grid.
  for (int i = 0; i < pyramid_.size(); i++) {
    const int num_grid_cells_per_dimension = 1 << (1 + i);
    pyramid_[i].setZero(num_grid_cells_per_dimension,
                        num_grid_cells_per_dimension);
  }
}

// Add a point to the visibility pyramid.
void VisibilityPyramid::AddPoint(const Eigen::Vector2d& point) {
  // Determine the grid cell of the point in the highest-resolution level of the
  // pyramid.
  int grid_cell_x = theia::Clamp(
      static_cast<int>(max_cells_in_dimension_ * point.x() / width_),
      0,
      max_cells_in_dimension_ - 1);
  int grid_cell_y = theia::Clamp(
      static_cast<int>(max_cells_in_dimension_ * point.y() / height_),
      0,
      max_cells_in_dimension_ - 1);

  // Go through the pyramid from fine to coarse and add the observation to the
  // occupancy grid.
  for (int i = pyramid_.size() - 1; i >= 0; --i) {
    pyramid_[i](grid_cell_x, grid_cell_y) += 1;

    // The next coarsest level of the pyramid will have half the number of grid
    // cells so we can use a simple bitshift to get the next pyramid level's
    // grid cells.
    grid_cell_x = grid_cell_x >> 1;
    grid_cell_y = grid_cell_y >> 1;
  }
}

// The score is accumulated by counting the number of grid cells in each level
// of the pyramid. The score for each level is weighted by the number of grid
// cells in that level of the pyramid. This scheme favors good spatial
// distribution at high resolutions.
int VisibilityPyramid::ComputeScore() const {
  int score = 0;
  for (int i = 0; i < pyramid_.size(); i++) {
    score += pyramid_[i].count() * pyramid_[i].size();
  }
  return score;
}

}  // namespace theia
