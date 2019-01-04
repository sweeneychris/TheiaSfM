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

// This solver is taken from http://cmp.felk.cvut.cz/~hellej1/ and implements
// the (non-minimal) six point
// (two-sided) radial distortion homography solver (H6_l1l2).
// Compared to the minimal five point solver (H5_l1l2) it yields only 2
// solutions and is almost twice as fast.

// This file was created by Steffen Urban (urbste@googlemail.com) or
// company address (steffen.urban@zeiss.com)
// January 2019

#ifndef THEIA_SFM_POSE_SIX_POINT_RADIAL_DISTORTION_HOMOGRAPHY_H_
#define THEIA_SFM_POSE_SIX_POINT_RADIAL_DISTORTION_HOMOGRAPHY_H_

#include <Eigen/Core>
#include <vector>

#include "theia/sfm/pose/six_point_radial_distortion_homography.h"

namespace theia {

struct RadialHomographyResult {
  Eigen::Matrix3d H;
  double l1;
  double l2;
};

// Input:
//   normalized_feature_points_left  - six normalized image positions  (inv(K)*p)
//   normalized_feature_points_right - six normalized image positions (inv(K)*p)
//   lmin - minimum radial distortion (can be used to speed up ransac loops)
//   lmax - maximum radial distortion (can be used to speed up ransac loops)
// Output: bool - returns true if at least one solution was found
// Return: RadialHomographyResult - struct that contains homography and a radial
// distortion parameter for each image (i.e. two-sided, called: H6_l1l2 in the paper)
bool SixPointRadialDistortionHomography(
    const std::vector<Eigen::Vector2d>& normalized_feature_points_left,
    const std::vector<Eigen::Vector2d>& normalized_feature_points_right,
    std::vector<RadialHomographyResult>* results, double lmin = -5.0,
    double lmax = 0.0);
}

#endif
