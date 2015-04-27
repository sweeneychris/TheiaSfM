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

#ifndef THEIA_SFM_BUNDLE_ADJUSTMENT_OPTIMIZE_RELATIVE_POSITION_WITH_KNOWN_ROTATION_H_
#define THEIA_SFM_BUNDLE_ADJUSTMENT_OPTIMIZE_RELATIVE_POSITION_WITH_KNOWN_ROTATION_H_

#include <Eigen/Core>
#include <vector>

#include "theia/matching/feature_correspondence.h"

namespace theia {

// Using known relative rotations, optimize the relative position that minimizes
// the epipolar constraint x2' * [t]_x * R * x1 = 0 for all
// correspondences. NOTE: the position is -R' * t and the rotations correspond
// to the absolute orientations of cameras 1 and 2.
//
// This algorithm is based on the robust pairwise translations estimation
// algorithm from "Robust Camera Location Estimation by Convex Programming" by
// Onur Ozyesil and Amit Singer (CVPR 2015).
bool OptimizeRelativePositionWithKnownRotation(
    const std::vector<FeatureCorrespondence>& correspondences,
    const Eigen::Vector3d& rotation1,
    const Eigen::Vector3d& rotation2,
    Eigen::Vector3d* relative_position);

}  // namespace theia

#endif  // THEIA_SFM_BUNDLE_ADJUSTMENT_OPTIMIZE_RELATIVE_POSITION_WITH_KNOWN_ROTATION_H_
