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

#ifndef THEIA_SFM_POSE_ESSENTIAL_MATRIX_UTILS_H_
#define THEIA_SFM_POSE_ESSENTIAL_MATRIX_UTILS_H_

#include <Eigen/Core>
#include <vector>

#include "theia/sfm/types.h"

namespace theia {

struct FeatureCorrespondence;

// Decomposes the essential matrix into the rotation R and translation t such
// that E can be any of the four candidate solutions: [rotation1 | translation],
// [rotation1 | -translation], [rotation2 | translation], [rotation2 |
// -translation].

void DecomposeEssentialMatrix(const Eigen::Matrix3d& essential_matrix,
                              Eigen::Matrix3d* rotation1,
                              Eigen::Matrix3d* rotation2,
                              Eigen::Vector3d* translation);

// Create an essential matrix from two projection matrices of the form [R|t].
void EssentialMatrixFromTwoProjectionMatrices(
    const Matrix3x4d& pose1,
    const Matrix3x4d& pose2,
    Eigen::Matrix3d* essential_matrix);

// Chooses the best pose of the 4 possible poses that can be computed from the
// essential matrix. The best pose is chosen as the pose that triangulates the
// most points in front of both cameras and the number of triangulated points is
// returned.
int GetBestPoseFromEssentialMatrix(
    const Eigen::Matrix3d& essential_matrix,
    const std::vector<FeatureCorrespondence>& normalized_correspondences,
    Eigen::Matrix3d* rotation,
    Eigen::Vector3d* position);

}  // namespace theia

#endif  // THEIA_SFM_POSE_ESSENTIAL_MATRIX_UTILS_H_
