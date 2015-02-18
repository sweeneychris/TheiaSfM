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

#ifndef THEIA_SFM_POSE_ESTIMATE_POSITIONS_LINEAR_H_
#define THEIA_SFM_POSE_ESTIMATE_POSITIONS_LINEAR_H_

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <unordered_map>
#include <vector>

#include "theia/sfm/reconstruction.h"
#include "theia/sfm/types.h"
#include "theia/sfm/view_triplet.h"
#include "theia/util/util.h"

namespace theia {

struct LinearPositionEstimatorOptions {
  int num_threads = 1;
};

// Estimates the camera position of views given global orientations and view
// triplets. The constraints formed by each triplet are used to create a sparse
// linear system to solve for the positions. This implementation closely follows
// "A Global Linear Method for Camera Pose Registration" by Jiang et al, ICCV
// 2013. Please see the paper for more details on the mathematics.
class LinearPositionEstimator {
 public:
  LinearPositionEstimator(const LinearPositionEstimatorOptions& options,
                          const Reconstruction& reconstruction);

  // Estimate the positions given triplets and global orientation estimates. No
  // filtering is done on the triplets to assure the triplets are of good
  // quality, so it is likely that you will want to filter the triplets before
  // calling this method.
  bool EstimatePositions(
      const std::vector<ViewTriplet>& triplets,
      const std::unordered_map<ViewId, Eigen::Vector3d>& orientation,
      std::unordered_map<ViewId, Eigen::Vector3d>* positions);

 private:
  // Computes the baselines ratios of each view pair in each triplet. These are
  // needed to estimate the position error.
  void ComputeBaselineRatios(const std::vector<ViewTriplet>& triplets,
                             std::vector<Eigen::Vector3d>* baselines);

  // Computes the baseline ratios of a single triplet.
  void ComputeBaselineRatioForTriplet(const ViewTriplet& triplet,
                                      Eigen::Vector3d* baseline);

  std::vector<TrackId> FindCommonTracks(const ViewId view_ids[3]);

  void GetTriangulatedPointDepths(const TwoViewInfo& info,
                                  const Feature& feature1,
                                  const Feature& feature2,
                                  double* depth1,
                                  double* depth2);

  Feature GetNormalizedFeature(const View& view, const TrackId track_id);

  // Sets up the linear system with the constraints that each triplet adds.
  void CreateLinearSystem(
      const std::vector<ViewTriplet>& triplets,
      const std::unordered_map<ViewId, Eigen::Vector3d>& orientations,
      std::vector<Eigen::Vector3d>& baselines,
      Eigen::SparseMatrix<double>* lhs,
      Eigen::VectorXd* rhs);

  const LinearPositionEstimatorOptions options_;
  const Reconstruction& reconstruction_;

  std::unordered_map<ViewId, int> num_triplets_for_view_;
  std::unordered_map<int, ViewId> linear_system_index_;

  DISALLOW_COPY_AND_ASSIGN(LinearPositionEstimator);
};

}  // namespace

#endif  // THEIA_SFM_POSE_ESTIMATE_POSITIONS_LINEAR_H_
