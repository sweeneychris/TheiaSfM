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

#include "theia/sfm/pose/compute_triplet_baseline_ratios.h"

#include <ceres/rotation.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <glog/logging.h>
#include <vector>

#include "theia/sfm/feature.h"
#include "theia/sfm/triangulation/triangulation.h"
#include "theia/sfm/view_triplet.h"

namespace theia {
namespace {

using Eigen::Vector3d;

// Triangulate the point and return the depth of the point relative to each
// view. Returns true if the depth was recovered successfully and false if the
// point could not be triangulated.
bool GetTriangulatedPointDepths(const TwoViewInfo& info,
                                const Vector3d& feature1,
                                const Vector3d& feature2,
                                double* depth1,
                                double* depth2) {
  static const double kMinTriangulationAngle = 2.0;

  const std::vector<Vector3d> origins = {Vector3d::Zero(), info.position_2};

  Eigen::Matrix3d rotation2;
  ceres::AngleAxisToRotationMatrix(
      info.rotation_2.data(), ceres::ColumnMajorAdapter3x3(rotation2.data()));
  const std::vector<Vector3d> directions = {feature1,
                                            rotation2.transpose() * feature2};

  // Make sure the rays have are viewed from a sufficient angle, otherwise the
  // depth computation is unstable.
  if (!SufficientTriangulationAngle(directions, kMinTriangulationAngle)) {
    return false;
  }

  Eigen::Vector4d point;
  if (!TriangulateMidpoint(origins, directions, &point)) {
    return false;
  }

  // Compute depths.
  const Vector3d point3d = point.hnormalized();
  *depth1 = point3d.norm();
  *depth2 = (point3d - info.position_2).norm();
  return true;
}

}  // namespace

void ComputeTripletBaselineRatios(const ViewTriplet& triplet,
                                  const std::vector<Feature>& feature1,
                                  const std::vector<Feature>& feature2,
                                  const std::vector<Feature>& feature3,
                                  Eigen::Vector3d* baseline) {
  CHECK_NOTNULL(baseline)->setZero();
  CHECK_EQ(feature1.size(), feature2.size())
      << "The feature containers must be the same size when computing the "
         "triplet baseline ratios.";
  CHECK_EQ(feature1.size(), feature3.size())
      << "The feature containers must be the same size when computing the "
         "triplet baseline ratios.";

  Eigen::Vector4d point12, point13, point23;
  double depth1_12, depth2_12, depth1_13, depth3_13, depth2_23, depth3_23;

  std::vector<double> baseline2, baseline3;
  baseline2.reserve(feature2.size());
  baseline3.reserve(feature3.size());
  for (int i = 0; i < feature1.size(); i++) {
    const Vector3d normalized_feature1 = feature1[i].homogeneous().normalized();
    const Vector3d normalized_feature2 = feature2[i].homogeneous().normalized();
    const Vector3d normalized_feature3 = feature3[i].homogeneous().normalized();
    if (!GetTriangulatedPointDepths(triplet.info_one_two,
                                    normalized_feature1,
                                    normalized_feature2,
                                    &depth1_12, &depth2_12)) {
      continue;
    }

    // Compute triangulation from views 1, 3.
    if (!GetTriangulatedPointDepths(triplet.info_one_three,
                                    normalized_feature1,
                                    normalized_feature3,
                                    &depth1_13, &depth3_13)) {
      continue;
    }

    // Compute triangulation from views 2, 3.
    if (!GetTriangulatedPointDepths(triplet.info_two_three,
                                    normalized_feature2,
                                    normalized_feature3,
                                    &depth2_23, &depth3_23)) {
      continue;
    }

    baseline2.emplace_back(depth1_12 / depth1_13);
    baseline3.emplace_back(depth2_12 / depth2_23);
  }
  CHECK_GT(baseline2.size(), 0) << "Could not compute the triplet baseline "
                                   "ratios. An inusfficient number of "
                                   "well-constrained 3D points were observed.";

  // Take the median as the baseline ratios.
  const int mid_index = baseline2.size() / 2;
  std::nth_element(baseline2.begin(),
                   baseline2.begin() + mid_index,
                   baseline2.end());
  std::nth_element(baseline3.begin(),
                   baseline3.begin() + mid_index,
                   baseline3.end());
  *baseline = Vector3d(1.0, baseline2[mid_index], baseline3[mid_index]);
}

}  // namespace theia
