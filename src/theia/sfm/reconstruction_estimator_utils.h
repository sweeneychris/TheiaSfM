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

#ifndef THEIA_SFM_RECONSTRUCTION_ESTIMATOR_UTILS_H_
#define THEIA_SFM_RECONSTRUCTION_ESTIMATOR_UTILS_H_

#include <Eigen/Core>
#include <unordered_map>

#include "theia/solvers/sample_consensus_estimator.h"
#include "theia/sfm/bundle_adjustment/bundle_adjustment.h"
#include "theia/sfm/reconstruction.h"
#include "theia/sfm/reconstruction_estimator.h"
#include "theia/sfm/pose/estimate_positions_nonlinear.h"
#include "theia/sfm/twoview_info.h"
#include "theia/sfm/view_graph/view_graph.h"

namespace theia {

// Sets the bundle adjustment optiosn from the reconstruction estimator options.
BundleAdjustmentOptions SetBundleAdjustmentOptions(
    const ReconstructionEstimatorOptions& options, const int num_views);

// Sets the ransac parameters from the reconstruction estimator options.
// NOTE: This does not set the error threshold since that is application
// specific. The caller must set this threshold.
RansacParameters SetRansacParameters(
    const ReconstructionEstimatorOptions& options);

// Initializes the focal length of each view. If EXIF data is available then the
// focal length is set using the EXIF value, otherwise it is set to
// 1.2 * max(image_width, image_height). This value is shown to be a decent
// initialization in the VisualSfM software.
void InitializeFocalLengthsFromImageSize(Reconstruction* reconstruction);

// Initializes focal length from the median focal length among all edges in the
// view graph.
void InitializeFocalLengthsFromMedian(const ViewGraph& view_graph,
                                      Reconstruction* reconstruction);

// Collects the relative rotations for each view pair into a simple map.
std::unordered_map<ViewIdPair, Eigen::Vector3d>
    RelativeRotationsFromTwoViewInfos(
        const std::unordered_map<ViewIdPair, TwoViewInfo>& view_pairs);


// Each view that has a rotation and position estimated has the Camera pose set.
void SetReconstructionFromEstimatedPoses(
    const std::unordered_map<ViewId, Eigen::Vector3d>& orientations,
    const std::unordered_map<ViewId, Eigen::Vector3d>& positions,
    Reconstruction* reconstruction);

// Outputs the ViewId of all estimated views in the reconstruction.
void GetEstimatedViewsFromReconstruction(const Reconstruction& reconstruction,
                                std::unordered_set<ViewId>* views);

// Outputs the TrackId of all estimated tracks in the reconstruction.
void GetEstimatedTracksFromReconstruction(const Reconstruction& reconstruction,
                                 std::unordered_set<TrackId>* tracks);

// Refine the relative translation estimates between view pairs by optimizing
// the epipolar constraint given the known rotation estimation.
void RefineRelativeTranslationsWithKnownRotations(
    const Reconstruction& reconstruction,
    const std::unordered_map<ViewId, Eigen::Vector3d>& orientations,
    std::unordered_map<ViewIdPair, TwoViewInfo>* view_pairs);

// Removes all features that have a reprojection error larger than the
// reprojection error threshold. Returns the number of features removed.
int RemoveOutlierFeatures(const double max_inlier_reprojection_error,
                          Reconstruction* reconstruction);

// Sets all tracks that are not seen by enough estimated views to unestimated.
// Returns the number of tracks set to unestimated.
int SetUnderconstrainedTracksToUnestimated(Reconstruction* reconstruction);

// Sets all vies that are not seen by enough estimated tracks to unestimated.
// Returns the number of views set to unestimated.
int SetUnderconstrainedViewsToUnestimated(Reconstruction* reconstruction);

}  // namespace theia

#endif  // THEIA_SFM_RECONSTRUCTION_ESTIMATOR_UTILS_H_
