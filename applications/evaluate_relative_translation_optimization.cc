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
//
// This program will evaluate the change in relative translation error produced
// by the relative translation optimization. The angular relative translation
// errors (compared the to ground truth relative translation) are shown before
// the optimization and after the optimization. The change in error is also
// output (negative and lower values mean the error decreased).xs

#include <Eigen/Core>
#include <glog/logging.h>
#include <gflags/gflags.h>
#include <theia/theia.h>

#include <memory>
#include <string>

DEFINE_string(
    matches, "",
    "Matches file that has been written with WriteMatchesAndGeometry");

DEFINE_string(reconstruction, "", "Reconstruction to use as ground truth.");

using theia::Reconstruction;
using theia::TrackId;
using theia::ViewId;
using theia::ViewIdPair;

// Accumulate all two view feature matches between the input views. The features
// are normalized according to the camera intrinsics.
void GetFeatureCorrespondences(
    const theia::View& view1,
    const theia::View& view2,
    std::vector<theia::FeatureCorrespondence>* matches) {
  Eigen::Matrix3d calibration1, calibration2;
  view1.Camera().GetCalibrationMatrix(&calibration1);
  view2.Camera().GetCalibrationMatrix(&calibration2);
  Eigen::Matrix3d inv_calibration1, inv_calibration2;
  bool view1_invertible, view2_invertible;
  double determinant;
  calibration1.computeInverseAndDetWithCheck(inv_calibration1, determinant,
                                             view1_invertible);
  calibration2.computeInverseAndDetWithCheck(inv_calibration2, determinant,
                                             view2_invertible);
  if (!view1_invertible || !view2_invertible) {
    LOG(FATAL) << "Calibration matrices are ill formed. Cannot optimize "
                  "epipolar constraints.";
    return;
  }

  const std::vector<theia::TrackId>& tracks = view1.TrackIds();
  for (const theia::TrackId track_id : tracks) {
    const theia::Feature* feature2 = view2.GetFeature(track_id);
    // If view 2 does not contain the current track then it cannot be a
    // correspondence.
    if (feature2 == nullptr) {
      continue;
    }

    theia::FeatureCorrespondence match;
    const theia::Feature* feature1 = view1.GetFeature(track_id);
    match.feature1 = *feature1;
    match.feature2 = *feature2;

    // Normalize for camera intrinsics.
    match.feature1 =
        (inv_calibration1 * match.feature1.homogeneous()).eval().hnormalized();
    match.feature2 =
        (inv_calibration2 * match.feature2.homogeneous()).eval().hnormalized();
    matches->emplace_back(match);
  }
}

double ComputeRelativeTranslationError(
    const Eigen::Vector3d& position1,
    const Eigen::Vector3d& position2,
    const Eigen::Matrix3d& rotation1,
    const Eigen::Vector3d& relative_translation) {
  const Eigen::Vector3d world_translation =
      rotation1 * (position2 - position1).normalized();
  return theia::RadToDeg(acos(
      theia::Clamp(relative_translation.dot(world_translation), -1.0, 1.0)));
}

void EvaluateTranslationOptimization(
    const std::string& matches_file,
    const theia::Reconstruction& reconstruction) {
  const std::vector<double> histogram_bins = {2,  5,   10,  15,  25,  50,
                                              90, 135, 180, 225, 270, 316};
  const std::vector<double> change_bins = {-50, -20, -10, -5, 0, 5, 10, 20, 50};
  theia::Histogram<double> before_hist(histogram_bins),
      after_hist(histogram_bins), change_hist(change_bins);

  std::vector<std::string> view_names;
  std::vector<theia::CameraIntrinsicsPrior> camera_intrinsics_prior;
  std::vector<theia::ImagePairMatch> matches;
  CHECK(theia::ReadMatchesAndGeometry(matches_file,
                                      &view_names,
                                      &camera_intrinsics_prior,
                                      &matches))
      << "Could not read matches from " << FLAGS_matches;

  // Collect relative translations.
  for (const auto& match : matches) {
    // Get the ViewIds from the names.
    const theia::ViewId view_id1 =
        reconstruction.ViewIdFromName(match.image1);
    const theia::ViewId view_id2 =
        reconstruction.ViewIdFromName(match.image2);
    const theia::View* view1 = reconstruction.View(view_id1);
    const theia::View* view2 = reconstruction.View(view_id2);
    if (view1 == nullptr || view2 == nullptr) {
      continue;
    }
    const theia::Camera& camera1 = view1->Camera();
    const theia::Camera& camera2 = view2->Camera();

    Eigen::Vector3d relative_position = match.twoview_info.position_2;
    const double translation_angular_error = ComputeRelativeTranslationError(
        camera1.GetPosition(),
        camera2.GetPosition(),
        camera1.GetOrientationAsRotationMatrix(),
        relative_position);
    before_hist.Add(translation_angular_error);

    // Optimize the relative position and recompute the error.
    std::vector<theia::FeatureCorrespondence> correspondences;
    GetFeatureCorrespondences(*view1, *view2, &correspondences);
    CHECK(theia::OptimizeRelativePositionWithKnownRotation(
        correspondences,
        camera1.GetOrientationAsAngleAxis(),
        camera2.GetOrientationAsAngleAxis(),
        &relative_position));

    const double optimized_translation_angular_error =
        ComputeRelativeTranslationError(
            camera1.GetPosition(),
            camera2.GetPosition(),
            camera1.GetOrientationAsRotationMatrix(),
            relative_position);
    after_hist.Add(optimized_translation_angular_error);

    change_hist.Add(optimized_translation_angular_error -
                    translation_angular_error);
  }
  LOG(INFO) << "Before histogram:\n" << before_hist.PrintString()
            << "\n\nAfter histogram:\n" << after_hist.PrintString()
            << "\n\nChange histogram:\n" << change_hist.PrintString();
}

int main(int argc, char* argv[]) {
  google::InitGoogleLogging(argv[0]);
  THEIA_GFLAGS_NAMESPACE::ParseCommandLineFlags(&argc, &argv, true);

  // Get the g.t. reconstruction first.
  std::unique_ptr<theia::Reconstruction> reconstruction(
      new theia::Reconstruction());
  CHECK(theia::ReadReconstruction(FLAGS_reconstruction, reconstruction.get()))
      << "Could not read reconstruction from " << FLAGS_reconstruction;

  LOG(INFO) << "Reading the matches.";
  EvaluateTranslationOptimization(FLAGS_matches, *reconstruction);

  return 0;
}
