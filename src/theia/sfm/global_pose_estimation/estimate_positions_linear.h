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

#ifndef THEIA_SFM_GLOBAL_POSE_ESTIMATION_ESTIMATE_POSITIONS_LINEAR_H_
#define THEIA_SFM_GLOBAL_POSE_ESTIMATION_ESTIMATE_POSITIONS_LINEAR_H_

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <unordered_map>
#include <vector>

#include "theia/sfm/reconstruction.h"
#include "theia/sfm/types.h"
#include "theia/sfm/view_triplet.h"
#include "theia/util/util.h"

namespace theia {

// Estimates the camera position of views given global orientations and view
// triplets. The constraints formed by each triplet are used to create a sparse
// linear system to solve for the positions. This implementation closely follows
// "A Global Linear Method for Camera Pose Registration" by Jiang et al, ICCV
// 2013. Please see the paper for more details on the mathematics.
class LinearPositionEstimator {
 public:
  struct Options {
    int num_threads = 1;

    // Maximum number of inverse power iterations to perform while extracting
    // the eigenvector corresponding to the smallest eigenvalue.
    int max_power_iterations = 1000;
    // The threshold at which to the iterative eigensolver method is considered
    // to be converged.
    double eigensolver_threshold = 1e-8;
  };

  LinearPositionEstimator(const Options& options,
                          const Reconstruction& reconstruction,
                          const std::vector<ViewTriplet>& triplets);

  // Estimate the positions given triplets and global orientation estimates. No
  // filtering is done on the triplets to assure the triplets are of good
  // quality, so it is likely that you will want to filter the triplets before
  // calling this method.
  bool EstimatePositions(
      const std::unordered_map<ViewId, Eigen::Vector3d>& orientation,
      std::unordered_map<ViewId, Eigen::Vector3d>* positions);

 private:
  // Returns the features as a unit-norm pixel ray after camera intrinsics
  // (i.e. focal length an principal point) have been removed.
  Feature GetNormalizedFeature(const View& view, const TrackId track_id);

  // Computes the relative baselines between three views in a triplet. The
  // baseline is estimated from the depths of triangulated 3D points.
  void ComputeBaselineRatioForTriplet(const ViewTriplet& triplet,
                                      Eigen::Vector3d* baseline);

  // Sets up the linear system with the constraints that each triplet adds.
  void CreateLinearSystem(
      const std::unordered_map<ViewId, Eigen::Vector3d>& orientations,
      const std::vector<Eigen::Vector3d>& baselines,
      Eigen::SparseMatrix<double>* constraint_matrix);

  // Positions are estimated from an eigenvector that is unit-norm with an
  // ambiguous sign. To ensure that the sign of the camera positions is correct,
  // we measure the relative translations from estimated camera positions and
  // compare that to the relative positions. If the sign is incorrect, we flip
  // the sign of all camera positions.
  void FlipSignOfPositionsIfNecessary(
      const std::unordered_map<ViewId, Eigen::Vector3d>& orientations,
      std::unordered_map<ViewId, Eigen::Vector3d>* positions);

  const Options options_;
  const Reconstruction& reconstruction_;
  const std::vector<ViewTriplet>& triplets_;

  // We keep one of the positions as constant to remove the ambiguity of the
  // origin of the linear system.
  static const int kConstantPositionIndex = -1;

  std::unordered_map<ViewId, int> num_triplets_for_view_;
  std::unordered_map<ViewId, int> linear_system_index_;

  DISALLOW_COPY_AND_ASSIGN(LinearPositionEstimator);
};

}  // namespace theia

#endif  // THEIA_SFM_GLOBAL_POSE_ESTIMATION_ESTIMATE_POSITIONS_LINEAR_H_
