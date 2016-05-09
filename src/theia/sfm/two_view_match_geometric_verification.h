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

#ifndef THEIA_SFM_TWO_VIEW_MATCH_GEOMETRIC_VERIFICATION_H_
#define THEIA_SFM_TWO_VIEW_MATCH_GEOMETRIC_VERIFICATION_H_

#include <vector>

#include "theia/matching/feature_correspondence.h"
#include "theia/matching/indexed_feature_match.h"
#include "theia/matching/keypoints_and_descriptors.h"
#include "theia/sfm/camera_intrinsics_prior.h"
#include "theia/sfm/camera/camera.h"
#include "theia/sfm/estimate_twoview_info.h"
#include "theia/sfm/twoview_info.h"
#include "theia/util/util.h"

namespace theia {

class TwoViewMatchGeometricVerification {
 public:
  struct Options {
    // Parameters for estimating the two view geometry.
    EstimateTwoViewInfoOptions estimate_twoview_info_options;

    // Minimum number of inlier matches in order to return true.
    int min_num_inlier_matches = 30;

    // TODO(csweeney): Guided matching.
    //
    // Perform guided matching to find more matches after initial geometry
    // estimation. Guided matching uses the current estimate for two-view
    // geometry to perform a constrained search along epipolar lines
    // corresponding to features. For a feature f in the first image, we search
    // for all features in the second image that lie near f's epipolar line
    // l. Among features close to l, we choose the feature with the smallest
    // descriptor distance as the match.
    // bool guided_matching = true;
    // For guided matching, features that are closer than this threshold to the
    // epipolar line will be considered for matching.
    // double guided_matching_max_distance_pixels = 2.0;

    // Bundle adjust the two view geometry using inliers.
    bool bundle_adjustment = true;

    // If performing bundle adjustment, the 3D points are only considered
    // inliers if the initial triangulation error is less than this. This value
    // is in pixels.
    double triangulation_max_reprojection_error = 15.0;

    // If performing bundle adjustment, the 3D points are only considered
    // inliers if the reprojection error after bundle adjustment is less than
    // this. This value is in pixels.
    double final_max_reprojection_error = 5.0;
  };

  TwoViewMatchGeometricVerification(
      const Options& options,
      const CameraIntrinsicsPrior& intrinsics1,
      const CameraIntrinsicsPrior& intrinsics2,
      const KeypointsAndDescriptors& features1,
      const KeypointsAndDescriptors& features2,
      const std::vector<IndexedFeatureMatch>& matches);

  // Perform 2-view geometric verification for the input. The verified matches
  // are returned along with the 2-view info. If the verification fails, false
  // is returned and the outputs are undefined.
  bool VerifyMatches(std::vector<FeatureCorrespondence>* verified_matches,
                     TwoViewInfo* twoview_info);

 private:
  // A helper method that creates a vector of FeatureCorrespondence from the
  // matches_ vector of match indices.
  void CreateCorrespondencesFromIndexedMatches(
      std::vector<FeatureCorrespondence>* correspondences);

  // Triangulates the current matches (given the cameras) and removes any
  // matches that do not have a small enough reprojection error after initial
  // triangulation (based on the options passed in).
  void TriangulatePoints(const Camera& camera1,
                         const Camera& camera2,
                         std::vector<Eigen::Vector4d>* triangulated_points);

  // Bundle adjusts the relative pose by triangulating 3D points. Points are
  // removed before and after bundle adjustment according to their reprojection
  // errors.
  bool BundleAdjustRelativePose(TwoViewInfo* twoview_info);

  // Estimates a homography and returns the number of inliers.
  int CountHomographyInliers();

  const Options options_;
  const CameraIntrinsicsPrior& intrinsics1_, intrinsics2_;
  const KeypointsAndDescriptors& features1_, features2_;

  // We keep a local copy of the matches so that we may add and remove matches
  // to it.
  std::vector<IndexedFeatureMatch> matches_;

  DISALLOW_COPY_AND_ASSIGN(TwoViewMatchGeometricVerification);
};

}  // namespace theia

#endif  // THEIA_SFM_TWO_VIEW_MATCH_GEOMETRIC_VERIFICATION_H_
