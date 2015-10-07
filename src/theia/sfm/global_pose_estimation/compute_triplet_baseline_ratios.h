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

#ifndef THEIA_SFM_GLOBAL_POSE_ESTIMATION_COMPUTE_TRIPLET_BASELINE_RATIOS_H_
#define THEIA_SFM_GLOBAL_POSE_ESTIMATION_COMPUTE_TRIPLET_BASELINE_RATIOS_H_

#include <Eigen/Core>
#include <vector>

#include "theia/sfm/feature.h"
#include "theia/sfm/view_triplet.h"

namespace theia {

// The baselines of the relative poses between views are estimated by
// triangulating features common to all 3 views. Based on the depth of the
// triangulated features, the baseline between views is recovered. The order of
// the features must be aligned such that feature1[i] corresponds to feature2[i]
// and feature3[i].
//
// The baselines returned in the 3-vector correspond to the baseline between
// views 1 and 2, between views 1 and 3, and between views 2 and 3 in that
// order.
//
// NOTE: The features must be normalized by the camera intrinsics (i.e.,
// principal point and focal length must be removed).
bool ComputeTripletBaselineRatios(const ViewTriplet& triplet,
                                  const std::vector<Feature>& feature1,
                                  const std::vector<Feature>& feature2,
                                  const std::vector<Feature>& feature3,
                                  Eigen::Vector3d* baseline);

}  // namespace theia

#endif  // THEIA_SFM_GLOBAL_POSE_ESTIMATION_COMPUTE_TRIPLET_BASELINE_RATIOS_H_
