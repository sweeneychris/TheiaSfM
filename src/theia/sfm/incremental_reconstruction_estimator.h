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

#ifndef THEIA_SFM_INCREMENTAL_RECONSTRUCTION_ESTIMATOR_H_
#define THEIA_SFM_INCREMENTAL_RECONSTRUCTION_ESTIMATOR_H_

#include <vector>
#include <unordered_map>

#include "theia/sfm/bundle_adjustment/bundle_adjustment.h"
#include "theia/sfm/estimate_track.h"
#include "theia/sfm/localize_view_to_reconstruction.h"
#include "theia/sfm/reconstruction_estimator.h"
#include "theia/sfm/reconstruction_estimator_options.h"
#include "theia/sfm/twoview_info.h"
#include "theia/sfm/types.h"
#include "theia/solvers/sample_consensus_estimator.h"
#include "theia/util/timer.h"
#include "theia/util/util.h"

namespace theia {

class Reconstruction;
class ViewGraph;

// Estimates the camera position and 3D structure of the scene using an
// incremental Structure from Motion approach. The method begins by first
// estimating the 3D structure and camera poses of 2 cameras based on their
// relative pose. Then additional cameras are added on sequentially and new 3D
// structure is estimated as new parts of the scene are observed. Bundle
// adjustment is repeatedly performed as more cameras are added to ensure high
// quality reconstructions and to avoid drift.
//
// The incremental SfM pipeline is as follows:
//   1) Choose an initial camera pair to reconstruct.
//   2) Estimate 3D structure of the scene.
//   3) Bundle adjustment on the 2-view reconstruction.
//   4) Localize a new camera to the current 3D points. Choose the camera that
//      observes the most 3D points currently in the scene.
//   5) Estimate new 3D structure.
//   6) Bundle adjustment if the model has grown by more than 5% since the last
//      bundle adjustment.
//   7) Repeat steps 4-6 until all cameras have been added.
//
// Incremental SfM is generally considered to be more robust than global SfM
// methods; hwoever, it requires many more instances of bundle adjustment (which
// is very costly) and so incremental SfM is not as efficient or scalable.
class IncrementalReconstructionEstimator : public ReconstructionEstimator {
 public:
  IncrementalReconstructionEstimator(
      const ReconstructionEstimatorOptions& options);

  ReconstructionEstimatorSummary Estimate(ViewGraph* view_graph,
                                          Reconstruction* reconstruction);

 private:
  // Choose two cameras to use as the seed for incremental reconstruction. These
  // cameras should observe 3D points that are well-conditioned. We determine
  // the conditioning of 3D points by examining the median viewing angle of the
  // correspondences between the views.
  bool ChooseInitialViewPair();

  // The best initial view pairs for incremental SfM are the ones that have a
  // lot of matches but sufficient baseline between them. One measurement for a
  // well-constrained baseline between cameras is the number of inliers when
  // estimating a homography (more inliers means it is less constrained i.e.,
  // bad). This method chooses the view pairs with more than
  // min_num_verified_matches that have the fewest homography inliers.
  void OrderViewPairsByInitializationCriterion(
      const int min_num_verified_matches,
      std::vector<ViewIdPair>* view_id_pairs);

  // Initialize the views based on the TwoViewInfo of the view pairs and set the
  // views as estimated.
  void InitializeCamerasFromTwoViewInfo(const ViewIdPair& view_ids);

  // Estimates all possible 3D points in the view. This is useful during
  // incremental SfM because we only need to triangulate points that were added
  // with new views.
  void EstimateStructure(const ViewId view_id);

  // The current percentage of cameras that have not been optimized by full BA.
  double UnoptimizedGrowthPercentage();

  // Performs partial bundle adjustment on the model. Only the k most recent
  // cameras (and the tracks observed in those views) are optimized.
  bool PartialBundleAdjustment();

  // Performs full bundle adjustment on the model.
  bool FullBundleAdjustment();

  // Chooses the next cameras to be localized according to which camera observes
  // the highest number of 3D points in the scene. This view is then localized
  // to using the calibrated or uncalibrated absolute pose algorithm.
  void FindViewsToLocalize(std::vector<ViewId>* views_to_localize);

  // Remove any features that have too high of reprojection errors or are not
  // well-constrained. Only the input features are checked for outliers.
  void RemoveOutlierTracks(const std::unordered_set<TrackId>& tracks_to_check,
                           const double max_reprojection_error_in_pixels);

  // Set any views that do not observe enough 3D points to unestimated, and
  // similarly set and tracks that are not observed by enough views to
  // unestimated.
  void SetUnderconstrainedAsUnestimated();

  ViewGraph* view_graph_;
  Reconstruction* reconstruction_;

  ReconstructionEstimatorOptions options_;
  BundleAdjustmentOptions bundle_adjustment_options_;
  RansacParameters ransac_params_;
  TrackEstimator::Options triangulation_options_;
  LocalizeViewToReconstructionOptions localization_options_;

  ReconstructionEstimatorSummary summary_;

  // A container to keep track of which views need to be localized.
  std::unordered_set<ViewId> unlocalized_views_;

  // An *ordered* container to keep track of which views have been added to the
  // reconstruction. This is used to determine which views are optimized during
  // partial BA.
  std::vector<ViewId> reconstructed_views_;

  // Indicates the number of views that have been optimized with full BA.
  int num_optimized_views_;

  DISALLOW_COPY_AND_ASSIGN(IncrementalReconstructionEstimator);
};

}  // namespace theia

#endif  // THEIA_SFM_INCREMENTAL_RECONSTRUCTION_ESTIMATOR_H_
