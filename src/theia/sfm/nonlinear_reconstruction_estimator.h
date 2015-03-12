// Copyright (C) 2014 The Regents of the University of California (Regents).
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

#ifndef THEIA_SFM_NONLINEAR_RECONSTRUCTION_ESTIMATOR_H_
#define THEIA_SFM_NONLINEAR_RECONSTRUCTION_ESTIMATOR_H_

#include "theia/sfm/bundle_adjustment/bundle_adjustment.h"
#include "theia/sfm/filter_view_pairs_from_relative_translation.h"
#include "theia/sfm/pose/estimate_positions_nonlinear.h"
#include "theia/sfm/reconstruction_estimator.h"
#include "theia/sfm/reconstruction_estimator_options.h"
#include "theia/sfm/twoview_info.h"
#include "theia/sfm/types.h"
#include "theia/solvers/sample_consensus_estimator.h"
#include "theia/util/util.h"

namespace theia {

class Reconstruction;
class ViewGraph;

// Estimates the camera position and 3D structure of the scene using nonlinear
// methods to estimate camera poses. First, rotation is estimated nonlinearlly
// then the position is estimated using a nonlinear optimization.
//
// The pipeline for estimating camera poses and structure is as follows:
//   1) Filter potentially bad pairwise geometries by enforcing a loop
//      constaint on rotations that form a triplet.
//   2) Initialize focal lengths.
//   3) Estimate the global rotation for each camera.
//   4) Remove any pairwise geometries where the relative rotation is not
//      consistent with the global rotation.
//   5) Optimize the relative translation given the known rotations.
//   6) Filter potentially bad relative translations with 1D SfM.
//   7) Estimate positions.
//   8) Estimate structure.
//   9) Bundle adjustment.
//   10) TODO: localize any unestimated cameras, retriangulate, and bundle
//      adjust.
//
// After each filtering step we remove any views which are no longer connected
// to the largest connected component in the view graph.
class NonlinearReconstructionEstimator : public ReconstructionEstimator {
 public:
  NonlinearReconstructionEstimator(
      const ReconstructionEstimatorOptions& options);

  ReconstructionEstimatorSummary Estimate(ViewGraph* view_graph,
                                          Reconstruction* reconstruction);

 private:
  bool FilterInitialViewGraph();
  void CalibrateCameras();
  void EstimateGlobalRotations();
  void FilterRotations();
  void OptimizePairwiseTranslations();
  void FilterRelativeTranslation();
  void EstimatePosition();
  void EstimateStructure();
  void BundleAdjustment();

  ViewGraph* view_graph_;
  Reconstruction* reconstruction_;

  ReconstructionEstimatorOptions options_;
  FilterViewPairsFromRelativeTranslationOptions translation_filter_options_;
  NonlinearPositionEstimatorOptions position_estimator_options_;
  BundleAdjustmentOptions bundle_adjustment_options_;
  RansacParameters ransac_params_;

  std::unordered_map<ViewIdPair, TwoViewInfo> view_pairs_;
  std::unordered_map<ViewId, Eigen::Vector3d> orientations_;
  std::unordered_map<ViewId, Eigen::Vector3d> positions_;

  DISALLOW_COPY_AND_ASSIGN(NonlinearReconstructionEstimator);
};

}  // namespace theia

#endif  // THEIA_SFM_NONLINEAR_RECONSTRUCTION_ESTIMATOR_H_
