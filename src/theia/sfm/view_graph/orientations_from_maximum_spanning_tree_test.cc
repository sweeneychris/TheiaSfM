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
#include <unordered_map>
#include <vector>

#include "gtest/gtest.h"
#include "theia/util/map_util.h"
#include "theia/sfm/view_graph/orientations_from_maximum_spanning_tree.h"
#include "theia/sfm/view_graph/view_graph.h"

namespace theia {

using Eigen::Matrix3d;
using Eigen::Vector3d;

void CreateViewsWithRandomOrientations(
    const int num_views,
    std::unordered_map<ViewId, Vector3d>* orientations) {
  (*orientations)[0] = Vector3d::Zero();
  for (int i = 1; i < num_views; i++) {
    (*orientations)[i] = Vector3d::Random();
  }
}

void RelativeRotationFromOrientations(const Vector3d& rotation1,
                                      const Vector3d& rotation2,
                                      Vector3d* relative_rotation) {
  Matrix3d orientation1, orientation2;
  ceres::AngleAxisToRotationMatrix(rotation1.data(), orientation1.data());
  ceres::AngleAxisToRotationMatrix(rotation2.data(), orientation2.data());
  const Matrix3d relative_rotation_mat =
      orientation2 * orientation1.transpose();
  ceres::RotationMatrixToAngleAxis(relative_rotation_mat.data(),
                                   relative_rotation->data());
}

TwoViewInfo CreateTwoViewInfo(
    const std::unordered_map<ViewId, Vector3d>& orientations,
    const ViewId view_id1,
    const ViewId view_id2) {
  TwoViewInfo info;
  RelativeRotationFromOrientations(FindOrDie(orientations, view_id1),
                                   FindOrDie(orientations, view_id2),
                                   &info.rotation_2);
  if (view_id1 > view_id2) {
    info.rotation_2 *= -1.0;
  }
  info.num_verified_matches = 100;
  return info;
}

void CreateViewGraph(
    const int num_edges,
    const std::unordered_map<ViewId, Vector3d>& orientations,
    ViewGraph* view_graph) {
  CHECK_GE(num_edges, orientations.size() - 1);

  // Add a skeletal graph.
  std::vector<ViewId> view_ids;
  view_ids.push_back(0);
  for (int i = 1; i < orientations.size(); i++) {
    const TwoViewInfo info = CreateTwoViewInfo(orientations, i - 1, i);
    view_graph->AddEdge(i - 1, i, info);
    view_ids.push_back(i);
  }

  // Add extra edges.
  while (view_graph->NumEdges() < num_edges) {
    std::random_shuffle(view_ids.begin(), view_ids.end());
    if (view_graph->GetEdge(view_ids[0], view_ids[1]) != nullptr) {
      continue;
    }
    const TwoViewInfo info =
        CreateTwoViewInfo(orientations, view_ids[0], view_ids[1]);
    view_graph->AddEdge(view_ids[0], view_ids[1], info);
  }
}

// Check that all edges in the view graph have the same relative rotations in
// the ground truth and estimated rotations.
void VerifyOrientations(
    const ViewGraph& view_graph,
    const std::unordered_map<ViewId, Vector3d>& gt_orientations,
    const std::unordered_map<ViewId, Vector3d>& orientations) {
  static const double kTolerance = 1e-12;
  for (const auto& edge : view_graph.GetAllEdges()) {
    const Vector3d& gt_orientation1 =
        FindOrDie(gt_orientations, edge.first.first);
    const Vector3d& gt_orientation2 =
        FindOrDie(gt_orientations, edge.first.second);
    Vector3d gt_relative_rotation;
    RelativeRotationFromOrientations(gt_orientation1,
                                     gt_orientation2,
                                     &gt_relative_rotation);

    const Vector3d& orientation1 = FindOrDie(orientations, edge.first.first);
    const Vector3d& orientation2 = FindOrDie(orientations, edge.first.second);
    Vector3d relative_rotation;
    RelativeRotationFromOrientations(orientation1,
                                     orientation2,
                                     &relative_rotation);

    EXPECT_LT((gt_relative_rotation - relative_rotation).norm(), kTolerance);
  }
}

void TestOrientationsFromViewGraph(const int num_views,
                                   const int num_edges) {
  std::unordered_map<ViewId, Vector3d> orientations;
  CreateViewsWithRandomOrientations(num_views, &orientations);
  ViewGraph view_graph;
  CreateViewGraph(num_edges, orientations, &view_graph);

  std::unordered_map<ViewId, Vector3d> estimated_orientations;
  OrientationsFromMaximumSpanningTree(view_graph, &estimated_orientations);
  VerifyOrientations(view_graph, orientations, estimated_orientations);
}

TEST(OrientationsFromViewGraph, SmallTest) {
  const int kNumViews = 4;
  const int kNumEdges = 6;
  TestOrientationsFromViewGraph(kNumViews, kNumEdges);
}

TEST(OrientationsFromViewGraph, LargeTest) {
  const int kNumViews = 50;
  const int kNumEdges = 100;
  TestOrientationsFromViewGraph(kNumViews, kNumEdges);
}

}  // namespace theia
