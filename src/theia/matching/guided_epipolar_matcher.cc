// Copyright (C) 2016 The Regents of the University of California (Regents).
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

#include "theia/matching/guided_epipolar_matcher.h"

#include <Eigen/Core>
#include <algorithm>
#include <limits>
#include <vector>
#include <unordered_map>
#include <unordered_set>

#include "theia/matching/keypoints_and_descriptors.h"
#include "theia/matching/distance.h"
#include "theia/matching/indexed_feature_match.h"
#include "theia/sfm/camera/camera.h"
#include "theia/sfm/pose/fundamental_matrix_util.h"
#include "theia/util/map_util.h"
#include "theia/util/hash.h"

namespace theia {

// Creates the grid structure for the fast epipolar lookup. This method must
// be called before calling GetMatches();
bool GuidedEpipolarMatcher::Initialize(
    const std::vector<IndexedFeatureMatch>& matches) {
  // Create a fast lookup for determining if a feature has already been matched.
  for (const IndexedFeatureMatch match : matches) {
    matched_features1_.insert(match.feature1_ind);
    matched_features2_.insert(match.feature2_ind);
  }

  image_grids_.emplace_back(
      ImageGrid(options_.guided_matching_max_distance_pixels, 0, 0));
  image_grids_.emplace_back(
      ImageGrid(options_.guided_matching_max_distance_pixels, 0.5, 0));
  image_grids_.emplace_back(
      ImageGrid(options_.guided_matching_max_distance_pixels, 0, 0.5));
  image_grids_.emplace_back(
      ImageGrid(options_.guided_matching_max_distance_pixels, 0.5, 0.5));

  // These helper vectors will be used to find the bounding box of the
  // features. This will help constrain the search along epipolar lines later.
  std::vector<double> x, y;
  x.reserve(features2_.keypoints.size());
  y.reserve(features2_.keypoints.size());
  // For each grid, add all features from features2.
  for (int i = 0; i < features2_.keypoints.size(); i++) {
    if (ContainsKey(matched_features2_, i)) {
      continue;
    }
    for (int j = 0; j < image_grids_.size(); j++) {
      image_grids_[j].AddFeature(i,
                                 features2_.keypoints[i].x(),
                                 features2_.keypoints[i].y());
    }
    x.emplace_back(features2_.keypoints[i].x());
    y.emplace_back(features2_.keypoints[i].y());
  }

  // Set the bounding box of the features.
  const auto& minmax_x = std::minmax_element(x.begin(), x.end());
  const auto& minmax_y = std::minmax_element(y.begin(), y.end());
  top_left_.x() = *minmax_x.first;
  top_left_.y() = *minmax_y.first;
  bottom_right_.x() = *minmax_x.second;
  bottom_right_.y() = *minmax_y.second;

  return true;
}

bool GuidedEpipolarMatcher::GetMatches(
    std::vector<IndexedFeatureMatch>* matches) {
  static const int kMinNumMatchesFound = 5;
  const int num_input_matches = matches->size();

  Initialize(*matches);

  // Compute the fundamental matrix based on the relative pose. The fundamental
  // matrix will map points in image 1 to lines in image 2.
  Eigen::Matrix<double, 3, 4> projection_matrix1, projection_matrix2;
  camera1_.GetProjectionMatrix(&projection_matrix1);
  camera2_.GetProjectionMatrix(&projection_matrix2);
  Eigen::Matrix3d fundamental_matrix;
  FundamentalMatrixFromProjectionMatrices(projection_matrix2.data(),
                                          projection_matrix1.data(),
                                          fundamental_matrix.data());

  // Perform guided matching for each feature.
  const int num_steps =
      static_cast<int>((top_left_ - bottom_right_).norm() /
                       options_.guided_matching_max_distance_pixels);
  L2 l2_distance;
  for (int i = 0; i < features1_.keypoints.size(); i++) {
    if (ContainsKey(matched_features1_, i)) {
      continue;
    }

    const Eigen::Vector3d point1(features1_.keypoints[i].x(),
                                 features1_.keypoints[i].y(),
                                 1.0);
    // Compute the epipolar line.
    Eigen::Vector3d epipolar_line = fundamental_matrix * point1;
    // Normalize the homogeneous line.
    epipolar_line /= epipolar_line.head<2>().norm();

    // Sample the epipolar line equally between the points where it intersects
    // the image boundaries.
    std::unordered_set<int> candidate_keypoints;
    for (int j = 0; j < num_steps; j++) {
      const Eigen::Vector2d sample_point =
          (j * top_left_ + (num_steps - j) * bottom_right_) / num_steps;

      // Find the cell center among all grids that is closest.
      std::vector<int> new_keypoints;
      FindClosestCellAndKeypoints(sample_point, &new_keypoints);
      candidate_keypoints.insert(new_keypoints.begin(), new_keypoints.end());
    }

    if (candidate_keypoints.size() < kMinNumMatchesFound) {
      continue;
    }

    // Compute the descriptor distances.
    std::vector<std::pair<float, int> > distances;
    distances.reserve(candidate_keypoints.size());
    for (const int match_index : candidate_keypoints) {
      distances.emplace_back(l2_distance(features1_.descriptors[i],
                                         features2_.descriptors[match_index]),
                             match_index);
    }

    // If the top 2 distance pass lowes ratio test then add the match to the
    // output.
    std::partial_sort(distances.begin(),
                      distances.begin() + 2,
                      distances.end());
    if (distances[0].first < distances[1].first * options_.lowes_ratio) {
      IndexedFeatureMatch match;
      match.feature1_ind = i;
      match.feature2_ind = distances[0].second;
      match.distance = distances[0].first;
      matches->emplace_back(match);
    }
  }

  const int num_added_matches = matches->size() - num_input_matches;
  LOG_IF(INFO, num_added_matches > 0)
      << "Guided matching added " << num_added_matches << " features to "
      << num_input_matches << " existing matches out of "
      << (std::min(features1_.keypoints.size(), features2_.keypoints.size()))
      << " possible matches.";
  return true;
}

void GuidedEpipolarMatcher::FindClosestCellAndKeypoints(
    const Eigen::Vector2d& point, std::vector<int>* new_keypoints) {
  int min_grid = 0;
  Eigen::Vector2i min_grid_center;
  double min_dist = std::numeric_limits<double>::max();

  // Find the grid cell with the closest center.
  for (int i = 0; i < 4; i++) {
    Eigen::Vector2i grid_center;
    image_grids_[i].GetClosestGridCenter(point.x(), point.y(), &grid_center);
    const double dist = (grid_center.cast<double>() - point).squaredNorm();
    if (dist < min_dist) {
      min_grid = i;
      min_dist = dist;
      std::swap(min_grid_center, grid_center);
    }
  }

  // Get the features from the closest grid cell.
  image_grids_[min_grid].GetFeaturesFromCell(min_grid_center, new_keypoints);
}

void GuidedEpipolarMatcher::ImageGrid::AddFeature(const int feature_index,
                                                  const double x,
                                                  const double y) {
  Eigen::Vector2i grid_center;
  GetClosestGridCenter(x, y, &grid_center);
  cells_[grid_center].emplace_back(feature_index);
}

void GuidedEpipolarMatcher::ImageGrid::GetClosestGridCenter(
    const double x, const double y, Eigen::Vector2i* grid_center) {
  grid_center->x() = static_cast<int>(
      std::floor(x / 2.0 * cell_size_ - cell_offset_x_) * cell_size_ +
      2.0 * cell_size_);
  grid_center->y() = static_cast<int>(
      std::floor(y / 2.0 * cell_size_ - cell_offset_y_) * cell_size_ +
      2.0 * cell_size_);
}

void GuidedEpipolarMatcher::ImageGrid::GetFeaturesFromCell(
    const Eigen::Vector2i& cell_center, std::vector<int>* feature_indices) {
  *feature_indices = FindWithDefault(cells_, cell_center, std::vector<int>());
}

}  // namespace theia
