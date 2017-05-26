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

#ifndef THEIA_SFM_GLOBAL_POSE_ESTIMATION_LINEAR_ROTATION_ESTIMATOR_H_
#define THEIA_SFM_GLOBAL_POSE_ESTIMATION_LINEAR_ROTATION_ESTIMATOR_H_

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <unordered_map>

#include "theia/sfm/global_pose_estimation/rotation_estimator.h"
#include "theia/sfm/types.h"
#include "theia/util/hash.h"

namespace theia {
class TwoViewInfo;

// Computes the orientation of views in a global frame given pairwise relative
// rotations between the views. This is done with a linear approximation to
// rotation averaging.
//
// This linear solution follows the method in "Robust Rotation and Translation
// Estimation in Multiview Geometry" by Martinec and Pajdla (CVPR 2007).
class LinearRotationEstimator : public RotationEstimator {
 public:
  LinearRotationEstimator() {}

  // Estimates the global orientations of all views based on an initial
  // guess. Returns true on successful estimation and false otherwise.
  bool EstimateRotations(
      const std::unordered_map<ViewIdPair, TwoViewInfo>& view_pairs,
      std::unordered_map<ViewId, Eigen::Vector3d>* global_orientations);

  // An alternative interface is to instead add relative rotation constraints
  // one by one with AddRelativeRotationConstraint, then call the
  // EstimateRotations interface below. This allows the caller to add multiple
  // constraints for the same view id pair, which may lead to more accurate
  // rotation estimates. Please see the following reference for an example of
  // how to obtain multiple constraints for pairs of views:
  //
  //   "Parallel Structure from Motion from Local Increment to Global Averaging"
  //   by Zhu et al (Arxiv 2017). https://arxiv.org/abs/1702.08601
  void AddRelativeRotationConstraint(const ViewIdPair& view_id_pair,
                                     const Eigen::Vector3d& relative_rotation);

  // Given the relative rotation constraints added with
  // AddRelativeRotationConstraint, this method returns the robust estimation of
  // global camera orientations. Like the method above, this requires an initial
  // estimate of the global orientations.
  bool EstimateRotations(
      std::unordered_map<ViewId, Eigen::Vector3d>* global_orientations);

 private:
  // Lookup map to keep track of the global orientation estimates by view id.
  std::unordered_map<ViewId, int> view_id_map_;

  // The sparse matrix is built up as new constraints are added.
  std::vector<Eigen::Triplet<double> > constraint_entries_;
};

}  // namespace theia

#endif  // THEIA_SFM_GLOBAL_POSE_ESTIMATION_LINEAR_ROTATION_ESTIMATOR_H_
