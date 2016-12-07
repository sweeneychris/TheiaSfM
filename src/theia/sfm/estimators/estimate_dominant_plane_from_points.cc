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
// Author: Benjamin Nuernberger (bnuernberger@cs.ucsb.edu)

#include "theia/sfm/estimators/estimate_dominant_plane_from_points.h"

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <limits>
#include <memory>
#include <vector>

#include "theia/sfm/create_and_initialize_ransac_variant.h"
#include "theia/sfm/pose/util.h"
#include "theia/solvers/estimator.h"
#include "theia/solvers/sample_consensus_estimator.h"
#include "theia/util/util.h"

namespace theia {
namespace {

using Eigen::Vector3d;

// An estimator for computing a dominant plane from a set of 3D points.
class DominantPlaneEstimator : public Estimator<Vector3d, Plane> {
 public:
  DominantPlaneEstimator() {}

  // 3 non-collinear points are needed to determine a plane.
  double SampleSize() const { return 3; }

  // Estimates candidate dominant planes from three 3D points.
  bool EstimateModel(const std::vector<Vector3d>& points,
                     std::vector<Plane>* planes) const {
    // If the points are collinear, there are no possible solutions.
    static const double kTolerance = 1e-6;
    const Vector3d a = points[1] - points[0];
    const Vector3d b = points[2] - points[0];
    const Vector3d cross = a.cross(b);
    if (cross.squaredNorm() < kTolerance) {
      VLOG(3)
          << "The 3 world points are collinear! No solution for a plane "
             "exists.";
      return false;
    }
    Plane plane;
    plane.point = points[0];
    plane.unit_normal = cross.normalized();

    planes->emplace_back(plane);
    return true;
  }

  // The error for a point given a plane model is the point-to-plane distance.
  double Error(const Vector3d& point, const Plane& plane) const {
    return std::abs(plane.unit_normal.dot(point - plane.point));
  }

 private:
  DISALLOW_COPY_AND_ASSIGN(DominantPlaneEstimator);
};

}  // namespace

bool EstimateDominantPlaneFromPoints(
    const RansacParameters& ransac_params,
    const RansacType& ransac_type,
    const std::vector<Vector3d>& points,
    Plane* plane,
    RansacSummary* ransac_summary) {
  DominantPlaneEstimator dominant_plane_estimator;
  std::unique_ptr<SampleConsensusEstimator<DominantPlaneEstimator> > ransac =
      CreateAndInitializeRansacVariant(ransac_type,
                                       ransac_params,
                                       dominant_plane_estimator);
  // Estimate the dominant plane.
  return ransac->Estimate(points, plane, ransac_summary);
}

}  // namespace theia
