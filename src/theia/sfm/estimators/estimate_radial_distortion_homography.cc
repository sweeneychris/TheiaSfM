// Copyright (C) 2019 The Regents of the University of California (Regents).
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

// This file was created by Steffen Urban (urbste@googlemail.com) or
// company address (steffen.urban@zeiss.com)
// May 2019

#include "theia/sfm/estimators/estimate_radial_distortion_homography.h"

#include <vector>

#include "theia/util/util.h"

namespace theia {
using Vector3d = Eigen::Vector3d;
using Vector2d = Eigen::Vector2d;
using Matrix3d = Eigen::Matrix3d;

// Vector3d Cam2CamH(const Matrix3d& H,
//                  const Vector3d &X)
//{
//    Vector3d Y = H * X;
//    Y /= Y(2);

//    return Y;
//}

// double RadialHomographyError(const RadialHomographyResult& radial_homography,
//                             const RadialDistortionFeatureCorrespondence
//                             image_correspondence)
//{
//    // set estimated radial distortion parameters
//    image_correspondence.camera_left->SetRadialDistortion(radial_homography.l1);
//    image_correspondence.camera_right->SetRadialDistortion(radial_homography.l2);

//    // unproject image points
//    const Vector3d ray_left =
//    image_correspondence.camera_left->ImageToCameraCoordinates(image_correspondence.feature_left);
//    const Vector3d ray_right =
//    image_correspondence.camera_right->ImageToCameraCoordinates(image_correspondence.feature_right);

//    const Vector3d ray_right_in_left = Cam2CamH(radial_homography.H,
//    ray_right);
//    const Vector3d ray_left_in_right = Cam2CamH(radial_homography.H.inverse(),
//    ray_left);

//    const Feature pt_left_projected  =
//    image_correspondence.camera_left->CameraToImageCoordinates(ray_right_in_left);
//    const Feature pt_right_projected =
//    image_correspondence.camera_right->CameraToImageCoordinates(ray_left_in_right);

//    const Feature diff_left = image_correspondence.feature_left -
//    pt_left_projected;
//    const Feature diff_right = image_correspondence.feature_right -
//    pt_right_projected;

//    const double squared_sum_left = diff_left.dot(diff_left);
//    const double squared_sum_right = diff_right.dot(diff_right);

//    return 0.5 * (squared_sum_left + squared_sum_right);
//}

// An estimator for computing the homography matrix from 6 feature
// correspondences.
class RadialHomographyMatrixEstimator
    : public Estimator<RadialDistortionFeatureCorrespondence,
                       RadialHomographyResult> {
 public:
  RadialHomographyMatrixEstimator() {}

  // 5 correspondences are needed to determine an essential matrix.
  double SampleSize() const { return 6; }

  // Estimates candidate essential matrices from correspondences.
  bool EstimateModel(
      const std::vector<RadialDistortionFeatureCorrespondence>& correspondences,
      std::vector<RadialHomographyResult>* results) const {
    std::vector<Vector2d> image_points_left(correspondences.size());
    std::vector<Vector2d> image_points_right(correspondences.size());
    for (int i = 0; i < correspondences.size(); ++i) {
      image_points_left[i] = correspondences[i].normalized_feature_left;
      image_points_right[i] = correspondences[i].normalized_feature_right;
    }

    return SixPointRadialDistortionHomography(
        image_points_left, image_points_right, results,
        correspondences[0].min_radial_distortion,
        correspondences[0].max_radial_distortion);
  }

  // The error for a correspondence given a model
  double Error(const RadialDistortionFeatureCorrespondence& correspondence,
               const RadialHomographyResult& radial_homography_result) const {
    return CheckRadialSymmetricError(
        radial_homography_result, correspondence.feature_left,
        correspondence.feature_right, correspondence.focal_length_estimate_left,
        correspondence.focal_length_estimate_right);
  }

 private:
  DISALLOW_COPY_AND_ASSIGN(RadialHomographyMatrixEstimator);
};

bool EstimateRadialHomographyMatrix(
    const RansacParameters& ransac_params, const RansacType& ransac_type,
    const std::vector<RadialDistortionFeatureCorrespondence>&
        normalized_correspondences,
    RadialHomographyResult* result, RansacSummary* ransac_summary) {
  RadialHomographyMatrixEstimator radial_homography_matrix_estimator;

  std::unique_ptr<SampleConsensusEstimator<RadialHomographyMatrixEstimator> >
      ransac = CreateAndInitializeRansacVariant(
          ransac_type, ransac_params, radial_homography_matrix_estimator);

  // Estimate essential matrix.
  return ransac->Estimate(normalized_correspondences, result, ransac_summary);
}
}
