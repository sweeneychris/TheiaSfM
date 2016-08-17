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

#include <ceres/rotation.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <unordered_map>
#include <vector>

#include "gtest/gtest.h"
#include "theia/math/util.h"
#include "theia/util/hash.h"
#include "theia/util/map_util.h"
#include "theia/util/random.h"
#include "theia/sfm/filter_view_pairs_from_relative_translation.h"
#include "theia/sfm/types.h"
#include "theia/sfm/view_graph/view_graph.h"

namespace theia {

namespace {

using Eigen::Matrix3d;
using Eigen::Vector3d;

RandomNumberGenerator rng(57);

void CreateViewsWithRandomPoses(
    const int num_views,
    std::unordered_map<ViewId, Vector3d>* orientations,
    std::unordered_map<ViewId, Vector3d>* positions) {
  (*orientations)[0] = Vector3d::Zero();
  (*positions)[0] = Vector3d::Zero();
  for (int i = 1; i < num_views; i++) {
    (*orientations)[i] = rng.RandVector3d();
    (*positions)[i] = rng.RandVector3d();
  }
}

TwoViewInfo CreateTwoViewInfo(
    const std::unordered_map<ViewId, Vector3d>& orientations,
    const std::unordered_map<ViewId, Vector3d>& positions,
    const ViewIdPair& view_id_pair) {
  TwoViewInfo info;

  Matrix3d orientation1, orientation2;
  ceres::AngleAxisToRotationMatrix(
      FindOrDie(orientations, view_id_pair.first).data(), orientation1.data());
  ceres::AngleAxisToRotationMatrix(
      FindOrDie(orientations, view_id_pair.second).data(), orientation2.data());
  const Matrix3d relative_rotation_mat =
      orientation2 * orientation1.transpose();
  ceres::RotationMatrixToAngleAxis(relative_rotation_mat.data(),
                                   info.rotation_2.data());

  const Vector3d position =
      (FindOrDie(positions, view_id_pair.second) -
       FindOrDie(positions, view_id_pair.first)).normalized();
  info.position_2 = orientation1 * position;

  return info;
}

void CreateValidViewPairs(
    const int num_valid_view_pairs,
    const std::unordered_map<ViewId, Vector3d>& orientations,
    const std::unordered_map<ViewId, Vector3d>& positions,
    ViewGraph* view_graph) {
  // Add a skeletal graph.
  std::vector<ViewId> view_ids;
  view_ids.push_back(0);
  for (int i = 1; i < orientations.size(); i++) {
    const ViewIdPair view_id_pair(i - 1, i);
    const TwoViewInfo info = CreateTwoViewInfo(orientations,
                                               positions,
                                               view_id_pair);
    view_graph->AddEdge(i - 1, i, info);
    view_ids.push_back(i);
  }

  // Add extra edges.
  while (view_graph->NumEdges() < num_valid_view_pairs) {
    std::random_shuffle(view_ids.begin(), view_ids.end());
    const ViewIdPair view_id_pair(view_ids[0], view_ids[1]);
    if (view_id_pair.first > view_id_pair.second ||
        view_graph->HasEdge(view_id_pair.first, view_id_pair.second)) {
      continue;
    }
    const TwoViewInfo info =
        CreateTwoViewInfo(orientations, positions, view_id_pair);
    view_graph->AddEdge(view_id_pair.first, view_id_pair.second, info);
  }
}

void CreateInvalidViewPairs(
    const int num_invalid_view_pairs,
    const std::unordered_map<ViewId, Vector3d>& orientations,
    const std::unordered_map<ViewId, Vector3d>& positions,
    ViewGraph* view_graph) {
  const int final_num_view_pairs =
      view_graph->NumEdges() + num_invalid_view_pairs;
  while (view_graph->NumEdges() < final_num_view_pairs) {
    // Choose a random view pair id.
    const ViewIdPair view_id_pair(rng.RandInt(0, orientations.size() - 1),
                                  rng.RandInt(0, orientations.size() - 1));
    if (view_id_pair.first >= view_id_pair.second ||
        view_graph->HasEdge(view_id_pair.first, view_id_pair.second)) {
      continue;
    }

    // Create a valid view pair.
    TwoViewInfo info =
        CreateTwoViewInfo(orientations, positions, view_id_pair);
    // Add a lot of noise to it.
    info.rotation_2 += Vector3d::Ones();
    info.position_2 = rng.RandVector3d().normalized();
    view_graph->AddEdge(view_id_pair.first, view_id_pair.second, info);
  }
}

void TestFilterViewPairsFromRelativeTranslation(
    const int num_views,
    const int num_valid_view_pairs,
    const int num_invalid_view_pairs) {
  std::unordered_map<ViewId, Vector3d> orientations;
  std::unordered_map<ViewId, Vector3d> positions;
  CreateViewsWithRandomPoses(num_views, &orientations, &positions);
  ViewGraph view_graph;
  CreateValidViewPairs(num_valid_view_pairs,
                       orientations,
                       positions,
                       &view_graph);
  CreateInvalidViewPairs(num_invalid_view_pairs,
                         orientations,
                         positions,
                         &view_graph);
  FilterViewPairsFromRelativeTranslationOptions options;
  options.rng = std::make_shared<RandomNumberGenerator>(169);
  FilterViewPairsFromRelativeTranslation(options, orientations, &view_graph);
  EXPECT_GE(view_graph.NumEdges(), num_valid_view_pairs);
}

}  // namespace

TEST(FilterViewPairsFromRelativeTranslation, LineTest) {
  static const int kNumViews = 4;
  static const int kValidViewPairs = 3;
  std::unordered_map<ViewId, Vector3d> orientations;
  std::unordered_map<ViewId, Vector3d> positions;
  for (int i = 0; i < kNumViews; i++) {
    orientations[i] = Vector3d::Zero();
    positions[i] = Vector3d(i, 0, 0);
  }

  ViewGraph view_graph;
  CreateValidViewPairs(kValidViewPairs,
                       orientations,
                       positions,
                       &view_graph);

  // Add two invalid view pairs.
  const ViewIdPair invalid_view_pair(0, 3);
  TwoViewInfo invalid_info;
  // Force the bad translations to be really bad.
  invalid_info.position_2 = Vector3d(-1, -1, -1).normalized();
  view_graph.AddEdge(0, 3, invalid_info);

  FilterViewPairsFromRelativeTranslationOptions options;
  options.rng = std::make_shared<RandomNumberGenerator>(199);
  options.translation_projection_tolerance = 0.1;
  FilterViewPairsFromRelativeTranslation(options, orientations, &view_graph);
  EXPECT_EQ(view_graph.NumEdges(), kValidViewPairs);
}

TEST(FilterViewPairsFromRelativeTranslation, NoBadRotations) {
  TestFilterViewPairsFromRelativeTranslation(10, 30, 0);
}

TEST(FilterViewPairsFromRelativeTranslation, FewBadRotations) {
  TestFilterViewPairsFromRelativeTranslation(10, 30, 5);
}

TEST(FilterViewPairsFromRelativeTranslation, ManyBadRotations) {
  TestFilterViewPairsFromRelativeTranslation(30, 100, 30);
}

}  // namespace theia
