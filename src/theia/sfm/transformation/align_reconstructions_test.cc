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

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <glog/logging.h>
#include <algorithm>
#include <utility>
#include <vector>
#include <string>

#include "gtest/gtest.h"

#include "theia/sfm/reconstruction.h"
#include "theia/sfm/types.h"
#include "theia/sfm/transformation/align_reconstructions.h"
#include "theia/sfm/transformation/transform_reconstruction.h"
#include "theia/util/stringprintf.h"

namespace theia {

Camera RandomCamera() {
  Camera camera;
  camera.SetPosition(10 * Eigen::Vector3d::Random());
  camera.SetOrientationFromAngleAxis(0.2 * Eigen::Vector3d::Random());
  camera.SetImageSize(1000, 1000);
  camera.SetFocalLength(800);
  camera.SetAspectRatio(1.0);
  camera.SetSkew(0.0);
  camera.SetPrincipalPoint(500.0, 500.0);
  return camera;
}

void BuildReconstructions(const int num_views,
                          const int num_tracks,
                          Reconstruction* reconstruction1,
                          Reconstruction* reconstruction2) {
  for (int i = 0; i < num_views; i++) {
    const std::string name = StringPrintf("%d", i);
    const ViewId view_id1 = reconstruction1->AddView(name);
    const ViewId view_id2 = reconstruction2->AddView(name);

    Camera camera = RandomCamera();
    *reconstruction1->MutableView(view_id1)->MutableCamera() = camera;
    *reconstruction2->MutableView(view_id2)->MutableCamera() = camera;
    reconstruction1->MutableView(view_id1)->SetEstimated(true);
    reconstruction2->MutableView(view_id2)->SetEstimated(true);
  }

  const std::vector<std::pair<ViewId, Feature> > track = {{0, Feature(0, 0)},
                                                          {1, Feature(0, 0)}};
  for (int i = 0; i < num_tracks; i++) {
    const TrackId track_id1 = reconstruction1->AddTrack(track);
    const TrackId track_id2 = reconstruction2->AddTrack(track);

    const Eigen::Vector4d point = Eigen::Vector4d::Random();
    *reconstruction1->MutableTrack(track_id1)->MutablePoint() = point;
    *reconstruction2->MutableTrack(track_id2)->MutablePoint() = point;
    reconstruction1->MutableTrack(track_id1)->SetEstimated(true);
    reconstruction2->MutableTrack(track_id2)->SetEstimated(true);
  }
}

void VerifyAlignment(const Reconstruction& reconstruction1,
                     const Reconstruction& reconstruction2) {
  static const double kTolerance = 1e-8;
  for (const ViewId view_id : reconstruction1.ViewIds()) {
    const Camera& camera1 = reconstruction1.View(view_id)->Camera();
    const Camera& camera2 = reconstruction2.View(view_id)->Camera();

    // Verify rotation.
    const Eigen::Vector3d rotation1 = camera1.GetOrientationAsAngleAxis();
    const Eigen::Vector3d rotation2 = camera2.GetOrientationAsAngleAxis();
    EXPECT_LT((rotation1 - rotation2).norm(), kTolerance);

    // Verify position.
    const Eigen::Vector3d position1 = camera1.GetPosition();
    const Eigen::Vector3d position2 = camera2.GetPosition();
    EXPECT_LT((position1 - position2).norm(), kTolerance)
        << "Position1 = " << position1.transpose()
        << " position2 = " << position2.transpose();
  }

  for (const TrackId track_id : reconstruction1.TrackIds()) {
    const Eigen::Vector3d point1 =
        reconstruction1.Track(track_id)->Point().hnormalized();
    const Eigen::Vector3d point2 =
        reconstruction2.Track(track_id)->Point().hnormalized();
    EXPECT_LT((point1 - point2).norm(), kTolerance)
        << "Point 1 = " << point1.transpose()
        << " point 2 = " << point2.transpose();
  }
}

void TestAlignReconstructions(const int num_views,
                              const int num_tracks,
                              const Eigen::Matrix3d rotation,
                              const Eigen::Vector3d translation,
                              const double scale) {
  Reconstruction reconstruction1, reconstruction2;
  BuildReconstructions(num_views,
                       num_tracks,
                       &reconstruction1,
                       &reconstruction2);
  TransformReconstruction(rotation, translation, scale, &reconstruction1);
  AlignReconstructions(reconstruction1, &reconstruction2);
  VerifyAlignment(reconstruction1, reconstruction2);
}

TEST(AlignReconstructions, Identity) {
  static const int kNumViews = 10;
  static const int kNumTracks = 20;
  const Eigen::Matrix3d rotation = Eigen::Matrix3d::Identity();
  const Eigen::Vector3d translation = Eigen::Vector3d::Zero();
  const double scale = 1.0;
  TestAlignReconstructions(kNumViews, kNumTracks, rotation, translation, scale);
}

TEST(AlignReconstructions, Scale) {
  static const int kNumViews = 10;
  static const int kNumTracks = 20;
  const Eigen::Matrix3d rotation = Eigen::Matrix3d::Identity();
  const Eigen::Vector3d translation = Eigen::Vector3d::Zero();
  const double scale = 5.3;
  TestAlignReconstructions(kNumViews, kNumTracks, rotation, translation, scale);
}

TEST(AlignReconstructions, Rotation) {
  static const int kNumViews = 10;
  static const int kNumTracks = 20;
  const Eigen::Vector3d rotation_aa = Eigen::Vector3d::Random();
  const Eigen::Matrix3d rotation =
      Eigen::AngleAxisd(rotation_aa.norm(), rotation_aa.normalized())
          .toRotationMatrix();
  const Eigen::Vector3d translation = Eigen::Vector3d::Zero();
  const double scale = 1.0;
  TestAlignReconstructions(kNumViews, kNumTracks, rotation, translation, scale);
}

TEST(AlignReconstructions, Translation) {
  static const int kNumViews = 10;
  static const int kNumTracks = 20;
  const Eigen::Matrix3d rotation = Eigen::Matrix3d::Identity();
  const Eigen::Vector3d translation = Eigen::Vector3d::Random();
  const double scale = 1.0;
  TestAlignReconstructions(kNumViews, kNumTracks, rotation, translation, scale);
}

TEST(AlignReconstructions, SimilarityTransformation) {
  static const int kNumViews = 10;
  static const int kNumTracks = 20;
  const Eigen::Vector3d rotation_aa = Eigen::Vector3d::Random();
  const Eigen::Matrix3d rotation =
      Eigen::AngleAxisd(rotation_aa.norm(), rotation_aa.normalized())
          .toRotationMatrix();
  const Eigen::Vector3d translation = Eigen::Vector3d::Random();
  const double scale = 4.2;
  TestAlignReconstructions(kNumViews, kNumTracks, rotation, translation, scale);
}


}  // namespace theia
