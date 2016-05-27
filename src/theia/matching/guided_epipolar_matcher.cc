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

#include "flann/flann.hpp"

#include "theia/matching/keypoints_and_descriptors.h"
#include "theia/matching/distance.h"
#include "theia/matching/indexed_feature_match.h"
#include "theia/sfm/camera/camera.h"
#include "theia/sfm/pose/fundamental_matrix_util.h"
#include "theia/util/map_util.h"
#include "theia/util/hash.h"
#include "theia/util/random.h"
#include "theia/util/timer.h"

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

  const double half_cell_width = options_.guided_matching_max_distance_pixels;
  const double offset = half_cell_width;
  image_grids_.emplace_back(ImageGrid(half_cell_width, 0, 0));
  image_grids_.emplace_back(ImageGrid(half_cell_width, offset, 0));
  image_grids_.emplace_back(ImageGrid(half_cell_width, 0, offset));
  image_grids_.emplace_back(ImageGrid(half_cell_width, offset, offset));

  // These helper vectors will be used to find the bounding box of the
  // features. This will help constrain the search along epipolar lines later.
  std::vector<double> x, y;
  x.reserve(features2_.keypoints.size());
  y.reserve(features2_.keypoints.size());
  // For each grid, add all features from features2.
  for (int i = 0; i < features2_.keypoints.size(); i++) {
    // TODO(csweeney): Test if the epipolar line of feature2 features is within
    // the image bounds of image 1. If not, we can skip the feature entirely.
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
  static const int kMinNumMatchesFound = 25;
  const int num_input_matches = matches->size();
  const double lowes_ratio_sq = options_.lowes_ratio * options_.lowes_ratio;

  Timer timer;
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

    // Find the points where the epipolar line intersects the bounding box of
    // the features in image 2. If the line is outside of the box (i.e. when we
    // cannot find two line endpoints), skip this feature for matching.
    //
    // TODO(csweeney): We should instead group features in image 1 by similar
    // epipolar lines (where line_endpoints[0] is within 2px) and use the same
    // KD tree for all features with this epipolar line. This should speed
    // things up significantly.
    std::vector<Eigen::Vector2d> line_endpoints;
    FindEpipolarLineIntersection(epipolar_line, &line_endpoints);
    if (line_endpoints.size() != 2) {
      continue;
    }
    const int num_steps =
        static_cast<int>((line_endpoints[1] - line_endpoints[0]).norm() /
                         options_.guided_matching_max_distance_pixels);

    // Sample the epipolar line equally between the points where it intersects
    // the features bounding box.
    std::vector<int> candidate_keypoint_indices;
    candidate_keypoint_indices.reserve(features2_.keypoints.size());
    for (int j = 0; j < num_steps; j++) {
      const Eigen::Vector2d sample_point =
          (j * line_endpoints[0] + (num_steps - j) * line_endpoints[1]) /
          num_steps;

      // Find the cell center among all grids that is closest.
      std::vector<int> new_keypoints;
      FindClosestCellAndKeypoints(sample_point, &new_keypoints);
      candidate_keypoint_indices.insert(candidate_keypoint_indices.end(),
                                        new_keypoints.begin(),
                                        new_keypoints.end());
    }

    // If we do not have enough features then the lowes ratio test below is not
    // meaningful. Add some random features here so that lowes ratio is more
    // informative of whether we have a good match or not.
    if (candidate_keypoint_indices.size() < kMinNumMatchesFound) {
      const int num_current_keypoints = candidate_keypoint_indices.size();
      for (int i = num_current_keypoints; i < kMinNumMatchesFound; i++) {
        candidate_keypoint_indices.emplace_back(
            RandInt(0, features2_.keypoints.size() - 1));
      }
    }

    // Remove duplicate entires in the candidate keypoints.
    std::sort(candidate_keypoint_indices.begin(),
              candidate_keypoint_indices.end());
    candidate_keypoint_indices.erase(
        std::unique(candidate_keypoint_indices.begin(),
                    candidate_keypoint_indices.end()),
        candidate_keypoint_indices.end());

    // Find the top 2 nearest neighbors.
    std::vector<std::pair<int, float> > nearest_neighbors;
    FindKNearestNeighbors(features1_.descriptors[i],
                          candidate_keypoint_indices,
                          &nearest_neighbors);

    // If the top 2 distance pass lowes ratio test then add the match to the
    // output.
    if (nearest_neighbors[0].second <
        nearest_neighbors[1].second * lowes_ratio_sq) {
      IndexedFeatureMatch match;
      match.feature1_ind = i;
      match.feature2_ind = nearest_neighbors[0].first;
      match.distance = nearest_neighbors[0].second;
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

void GuidedEpipolarMatcher::FindEpipolarLineIntersection(
    const Eigen::Vector3d& epipolar_line,
    std::vector<Eigen::Vector2d>* lines) {
  const double y_of_left_intersection =
      -(epipolar_line.z() + epipolar_line.x() * top_left_.x()) /
      epipolar_line.y();
  if (y_of_left_intersection >= top_left_.y() &&
      y_of_left_intersection <= bottom_right_.y()) {
    lines->emplace_back(Eigen::Vector2d(top_left_.x(), y_of_left_intersection));
  }

  const double x_of_top_intersection =
      -(epipolar_line.z() + epipolar_line.y() * top_left_.y()) /
      epipolar_line.x();
  if (x_of_top_intersection >= top_left_.x() &&
      x_of_top_intersection <= bottom_right_.x()) {
    lines->emplace_back(Eigen::Vector2d(x_of_top_intersection, top_left_.y()));
  }

  const double y_of_right_intersection =
      -(epipolar_line.z() + epipolar_line.x() * bottom_right_.x()) /
      epipolar_line.y();
  if (y_of_right_intersection >= top_left_.y() &&
      y_of_right_intersection <= bottom_right_.y()) {
    lines->emplace_back(
        Eigen::Vector2d(bottom_right_.x(), y_of_right_intersection));
  }

  const double x_of_bottom_intersection =
      -(epipolar_line.z() + epipolar_line.y() * bottom_right_.y()) /
      epipolar_line.x();
  if (x_of_bottom_intersection >= top_left_.x() &&
      x_of_bottom_intersection <= bottom_right_.x()) {
    lines->emplace_back(
        Eigen::Vector2d(x_of_bottom_intersection, bottom_right_.y()));
  }
}

void GuidedEpipolarMatcher::FindKNearestNeighbors(
    const Eigen::VectorXf& query_descriptor,
    const std::vector<int>& candidate_keypoint_indices,
    std::vector<std::pair<int, float> >* matches) {
  static const int kNumQueryDescriptors = 1;
  static const int kNumNearestNeighbors = 2;
  static const int kMinNumLeafsVisited = 50;

  // Gather the candidate matching descriptors.
  Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
      candidate_descriptors(candidate_keypoint_indices.size(),
                            query_descriptor.size());
  for (int i = 0; i < candidate_keypoint_indices.size(); i++) {
    const int match_index = candidate_keypoint_indices[i];
    candidate_descriptors.row(i) = features2_.descriptors[match_index];
  }

  // Create the searchable KD-tree with FLANN.
  flann::Matrix<float> flann_descriptors(candidate_descriptors.data(),
                                         candidate_descriptors.rows(),
                                         candidate_descriptors.cols());

  flann::Index<flann::L2<float> > flann_kd_tree(
      flann_descriptors, flann::KDTreeSingleIndexParams());
  flann_kd_tree.buildIndex();

  // Query the KD-tree to get the top 2 nearest neighbors.
  Eigen::RowVectorXf eig_query_descriptor = query_descriptor;
  const flann::Matrix<float> flann_query(eig_query_descriptor.data(),
                                         kNumQueryDescriptors,
                                         eig_query_descriptor.size());
  std::vector<std::vector<int> > nn_indices;
  std::vector<std::vector<float> > nn_distances;
  const int max_leafs_to_check =
      std::max(static_cast<int>(candidate_descriptors.rows() * 0.2),
               kMinNumLeafsVisited);
  flann_kd_tree.knnSearch(flann_query, nn_indices, nn_distances,
                          kNumNearestNeighbors,
                          flann::SearchParams(max_leafs_to_check));

  // Output the top 2 matches.
  matches->emplace_back(candidate_keypoint_indices[nn_indices[0][0]],
                        nn_distances[0][0]);
  matches->emplace_back(candidate_keypoint_indices[nn_indices[0][1]],
                        nn_distances[0][1]);
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

GuidedEpipolarMatcher::ImageGrid::ImageGrid(const double cell_size,
                                            const double cell_offset_x,
                                            const double cell_offset_y)
    : cell_size_(cell_size),
      cell_offset_x_(cell_offset_x),
      cell_offset_y_(cell_offset_y) {
  cells_.reserve(10000);
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
      std::floor((x - cell_offset_x_) / (2.0 * cell_size_)) * 2.0 * cell_size_ +
      cell_size_ + cell_offset_x_);
  grid_center->y() = static_cast<int>(
      std::floor((y - cell_offset_y_) / (2.0 * cell_size_)) * 2.0 * cell_size_ +
      + cell_size_ + cell_offset_y_);
}

void GuidedEpipolarMatcher::ImageGrid::GetFeaturesFromCell(
    const Eigen::Vector2i& cell_center, std::vector<int>* feature_indices) {
  *feature_indices = FindWithDefault(cells_, cell_center, std::vector<int>());
}

}  // namespace theia
