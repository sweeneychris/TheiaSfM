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

#ifndef THEIA_SFM_POSE_FIVE_POINT_RELATIVE_POSE_H_
#define THEIA_SFM_POSE_FIVE_POINT_RELATIVE_POSE_H_

#include <Eigen/Core>
#include <vector>

#include "theia/alignment/alignment.h"

namespace theia {

// Computes the relative pose between two cameras using 5 corresponding
// points. Algorithm is implemented based on "H. Stewénius, C. Engels, and
// D. Nistér. Recent developments on direct relative orientation". ISPRS Journal
// of Photogrammetry and Remote Sensing, 2006. The relative pose is computed
// such that y * E * x = 0, where E = t_x * R and t_x is the cross product
// matrix of t. This implementation is proven to be more stable than the fast 5
// pt. algorithm by Nister.
//
// Params:
//   image1_points: Location of features on the image plane of image 1.
//   image2_points: Location of features on the image plane of image 2.
// Return: true if a valid solution was found.
//
// NOTE: At least 5 points must be supplied, but a non-minimal estimate will be
// computed if more than five are supplied.
bool FivePointRelativePose(const std::vector<Eigen::Vector2d>& image1_points,
                           const std::vector<Eigen::Vector2d>& image2_points,
                           std::vector<Eigen::Matrix3d>* essential_matrices);
}  // namespace theia

#endif  // THEIA_SFM_POSE_FIVE_POINT_RELATIVE_POSE_H_
