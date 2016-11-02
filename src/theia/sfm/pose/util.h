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

#ifndef THEIA_SFM_POSE_UTIL_H_
#define THEIA_SFM_POSE_UTIL_H_

#include <Eigen/Core>
#include "theia/alignment/alignment.h"

namespace theia {

// Calculates Sampson distance for two correspondances and an essential or
// fundamental matrix by eq. 11.9 in Hartley and Zisserman. For an E or F
// that is defined such that y^t * E * x = 0
double SquaredSampsonDistance(const Eigen::Matrix3d& F,
                              const Eigen::Vector2d& x,
                              const Eigen::Vector2d& y);

// Returns the cross product matrix of a vector: if cross_vec = [x y z] then
//                        [ 0  -z   y]
// cross product matrix = [ z   0  -y]
//                        [-y   x   0]
Eigen::Matrix3d CrossProductMatrix(const Eigen::Vector3d& cross_vec);

// Given a 2xN matrix of image points (of the form [x, y]), this method
// calculates the matrix that will shift the points so that the centroid is at
// the origin and the average distance from the centroid is sqrt(2). Returns the
// transformation matrix and the transformed points.
bool NormalizeImagePoints(
    const std::vector<Eigen::Vector2d>& image_points,
    std::vector<Eigen::Vector2d>* normalized_image_points,
    Eigen::Matrix3d* normalization_matrix);

// Projects a 3x3 matrix to the rotation matrix in SO3 space with the closest
// Frobenius norm. For a matrix with an SVD decomposition M = USV, the nearest
// rotation matrix is R = UV'.
Eigen::Matrix3d ProjectToRotationMatrix(const Eigen::Matrix3d& matrix);

}  // namespace theia

#endif  // THEIA_SFM_POSE_UTIL_H_
