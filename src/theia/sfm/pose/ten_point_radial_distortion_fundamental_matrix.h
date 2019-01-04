// Copyright (C) 2019 The Regents of the University of California (Regents).
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

// This file was created by Steffen Urban (urbste@googlemail.com) or
// company address (steffen.urban@zeiss.com)
// January 2019

// This solver is taken from http://cmp.felk.cvut.cz/~hellej1/ and implements
// the ten point (two-sided) radial distortion fundamental matrix solver
// (H10_l1l2).
// The corresponding paper is:
// "Efficient Solution to the Epipolar Geometry for Radially Distorted Cameras",
// Zuzana Kukelova et al. ICCV 2015

#ifndef THEIA_SFM_POSE_TEN_POINT_RADIAL_DISTORTION_FUNDAMENTAL_MATRIX_H_
#define THEIA_SFM_POSE_TEN_POINT_RADIAL_DISTORTION_FUNDAMENTAL_MATRIX_H_

#include <Eigen/Core>
#include <vector>

#include "theia/sfm/pose/ten_point_radial_distortion_fundamental_matrix.h"

namespace theia {

struct RadialFundamentalMatrixResult {
  Eigen::Matrix3d F;
  double l1;
  double l2;
};

// Input:
//   normalized_feature_points_left  - ten normalized image positions (inv(K)*p)
//   normalized_feature_points_right - ten normalized image positions (inv(K)*p)
//   lmin - minimum radial distortion (can be used to speed up ransac loops)
//   lmax - maximum radial distortion (can be used to speed up ransac loops)
// Output: bool - returns true if at least one solution was found
// Return: RadialFundamentalMatrixResult - struct that contains fundamental
// matrix and a radial
// distortion parameter for each image (i.e. two-sided, called: H10_l1l2 in the
// paper)
bool TenPointRadialDistortionFundamentalMatrix(
    const std::vector<Eigen::Vector2d>& normalized_feature_points_left,
    const std::vector<Eigen::Vector2d>& normalized_feature_points_right,
    std::vector<RadialFundamentalMatrixResult>* results, double lmin = -5.0,
    double lmax = 0.0);
}

#endif
