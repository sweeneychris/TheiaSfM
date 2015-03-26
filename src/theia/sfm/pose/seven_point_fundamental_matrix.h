// Copyright (C) 2015 The Regents of the University of California (Regents).
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

#ifndef THEIA_SFM_POSE_SEVEN_POINT_FUNDAMENTAL_MATRIX_H_
#define THEIA_SFM_POSE_SEVEN_POINT_FUNDAMENTAL_MATRIX_H_

#include <Eigen/Core>
#include <vector>

namespace theia {
// Computes the Fundamental Matrix from seven points. This is a minimal solver
// for the the fundmanetal matrix and will return either 1 or 3 valid solutions
// (complex solutions are discarded).
//
// Params:
//   image1_points: 7 image points from one image.
//   image2_points: 7 image points from a second image.
//   fundamental_matrix: the estimated fundamental matrix such that
//     x2^t * F * x1 = 0 for points x1 in image1_points and x2 in
//     image2_points.
bool SevenPointFundamentalMatrix(
    const std::vector<Eigen::Vector2d>& image1_points,
    const std::vector<Eigen::Vector2d>& image2_points,
    std::vector<Eigen::Matrix3d>* fundamental_matrices);

}  // namespace theia

#endif  // THEIA_SFM_POSE_SEVEN_POINT_FUNDAMENTAL_MATRIX_H_
