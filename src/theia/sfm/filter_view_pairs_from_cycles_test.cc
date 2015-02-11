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

#include <ceres/rotation.h>
#include <Eigen/Core>
#include <unordered_map>
#include <vector>

#include "gtest/gtest.h"
#include "theia/util/hash.h"
#include "theia/util/map_util.h"
#include "theia/util/random.h"
#include "theia/sfm/filter_view_pairs_from_cycles.h"
#include "theia/sfm/types.h"

namespace theia {

namespace {

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
    const ViewIdPair& view_id_pair) {
  TwoViewInfo info;
  RelativeRotationFromOrientations(FindOrDie(orientations, view_id_pair.first),
                                   FindOrDie(orientations, view_id_pair.second),
                                   &info.rotation_2);
  return info;
}

void CreateValidViewPairs(
    const int num_valid_view_pairs,
    const std::unordered_map<ViewId, Vector3d>& orientations,
    std::unordered_map<ViewIdPair, TwoViewInfo>* view_pairs) {
  // Add a skeletal graph.
  std::vector<ViewId> view_ids;
  view_ids.push_back(0);
  for (int i = 1; i < orientations.size(); i++) {
    const ViewIdPair view_id_pair(i - 1, i);
    const TwoViewInfo info = CreateTwoViewInfo(orientations, view_id_pair);
    (*view_pairs)[view_id_pair] = info;
    view_ids.push_back(i);
  }

  // Add extra edges.
  while (view_pairs->size() < num_valid_view_pairs) {
    std::random_shuffle(view_ids.begin(), view_ids.end());
    const ViewIdPair view_id_pair(view_ids[0], view_ids[1]);
    if (view_id_pair.first > view_id_pair.second ||
        ContainsKey(*view_pairs, view_id_pair)) {
      continue;
    }
    const TwoViewInfo info =
        CreateTwoViewInfo(orientations, view_id_pair);
    (*view_pairs)[view_id_pair] = info;
  }
}

void CreateInvalidViewPairs(
    const int num_invalid_view_pairs,
    const std::unordered_map<ViewId, Vector3d>& orientations,
    std::unordered_map<ViewIdPair, TwoViewInfo>* view_pairs) {
  InitRandomGenerator();

  const int final_num_view_pairs = view_pairs->size() + num_invalid_view_pairs;
  while (view_pairs->size() < final_num_view_pairs) {
    // Choose a random view pair id.
    const ViewIdPair view_id_pair(RandInt(0, orientations.size() - 1),
                                  RandInt(0, orientations.size() - 1));
    if (view_id_pair.first == view_id_pair.second ||
        ContainsKey(*view_pairs, view_id_pair)) {
      continue;
    }

    // Create a valid view pair.
    (*view_pairs)[view_id_pair] = CreateTwoViewInfo(orientations, view_id_pair);
    // Add a lot of noise to it.
    (*view_pairs)[view_id_pair].rotation_2 += Vector3d::Ones();
  }
}

void TestFilterViewPairsFromCycles(const int num_views,
                                   const int num_valid_view_pairs,
                                   const int num_invalid_view_pairs) {
  static const double kMaxRelativeRotationDifferenceDegrees = 2.0;
  std::unordered_map<ViewId, Vector3d> orientations;
  CreateViewsWithRandomOrientations(num_views, &orientations);
  std::unordered_map<ViewIdPair, TwoViewInfo> view_pairs;
  CreateValidViewPairs(num_valid_view_pairs, orientations, &view_pairs);
  CreateInvalidViewPairs(num_invalid_view_pairs, orientations, &view_pairs);
  FilterViewPairsFromCycles(kMaxRelativeRotationDifferenceDegrees, &view_pairs);
  EXPECT_EQ(view_pairs.size(), num_valid_view_pairs);
}

}  // namespace

TEST(FilterViewPairsFromCycles, NoBadRotations) {
  TestFilterViewPairsFromCycles(10, 30, 0);
}

TEST(FilterViewPairsFromCycles, FewBadRotations) {
  TestFilterViewPairsFromCycles(10, 30, 5);
}

TEST(FilterViewPairsFromCycles, ManyBadRotations) {
  TestFilterViewPairsFromCycles(10, 30, 30);
}

}  // namespace theia
