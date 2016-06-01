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
#include <stdint.h>
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

namespace theia {
namespace {

// Encodes the line endpoints into an uint64_t for fast sorting.
uint64_t EncodeLineEndpoints(const std::vector<Eigen::Vector2d>& endpoints) {
  uint64_t encoded_endpoint = 0;
  const uint16_t x1 = static_cast<uint16_t>(endpoints[0].x());
  const uint16_t y1 = static_cast<uint16_t>(endpoints[0].y());
  const uint16_t x2 = static_cast<uint16_t>(endpoints[1].x());
  const uint16_t y2 = static_cast<uint16_t>(endpoints[1].y());
  encoded_endpoint |= static_cast<uint64_t>(x1) << 48;
  encoded_endpoint |= static_cast<uint64_t>(y1) << 32;
  encoded_endpoint |= static_cast<uint64_t>(x2) << 16;
  encoded_endpoint |= static_cast<uint64_t>(y2) << 0;
  return encoded_endpoint;
}

void DecodeLineEndpoints(const uint64_t encoded_endpoint,
                         std::vector<Eigen::Vector2d>* endpoints) {
  // Create a bitmask to isolate the last 16 bits.
  uint64_t s = 0;
  s = ~s;
  s = s >> 48;

  const uint16_t x1 = static_cast<uint16_t>(s & (encoded_endpoint >> 48));
  const uint16_t y1 = static_cast<uint16_t>(s & (encoded_endpoint >> 32));
  const uint16_t x2 = static_cast<uint16_t>(s & (encoded_endpoint >> 16));
  const uint16_t y2 = static_cast<uint16_t>(s & (encoded_endpoint >> 0));

  endpoints->emplace_back(static_cast<double>(x1), static_cast<double>(y1));
  endpoints->emplace_back(static_cast<double>(x2), static_cast<double>(y2));
}

}  // namespace

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
  const int num_input_matches = matches->size();
  const double lowes_ratio_sq = options_.lowes_ratio * options_.lowes_ratio;

  Initialize(*matches);

  // Group all epipolar lines.
  std::vector<EpilineGroup> epiline_groups;
  GroupEpipolarLines(&epiline_groups);

  for (const EpilineGroup& epiline_group : epiline_groups) {
    // Find all features close to this epiline.
    std::vector<int> candidate_keypoint_indices;
    FindFeaturesNearEpipolarLines(epiline_group, &candidate_keypoint_indices);

    // Build a KD-tree for those features and query the tree to get the nearest
    // neighbor matches.
    std::vector<std::vector<float> > nn_distances;
    std::vector<std::vector<int> > nn_indices;
    FindKNearestNeighbors(epiline_group.features, candidate_keypoint_indices,
                          &nn_distances, &nn_indices);

    // For each of the nearest neighbor matches, check the Lowes ratio to
    // determine if the match is valid.
    for (int i = 0; i < nn_distances.size(); i++) {
      // If the top 2 distance pass lowes ratio test then add the match to the
      // output.
      if (nn_distances[i][0] < nn_distances[i][1] * lowes_ratio_sq) {
        IndexedFeatureMatch match;
        match.feature1_ind = epiline_group.features[i];
        match.feature2_ind = nn_indices[i][0];
        match.distance = nn_distances[i][0];
        matches->emplace_back(match);
      }
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


void GuidedEpipolarMatcher::GroupEpipolarLines(
    std::vector<EpilineGroup>* epiline_groups) {
  // Compute the fundamental matrix based on the relative pose. The fundamental
  // matrix will map points in image 1 to lines in image 2.
  Eigen::Matrix<double, 3, 4> projection_matrix1, projection_matrix2;
  camera1_.GetProjectionMatrix(&projection_matrix1);
  camera2_.GetProjectionMatrix(&projection_matrix2);
  Eigen::Matrix3d fundamental_matrix;
  FundamentalMatrixFromProjectionMatrices(projection_matrix2.data(),
                                          projection_matrix1.data(),
                                          fundamental_matrix.data());

  // Sort the endpoints by encoding x1, y1, x2, y2 into a long long int. This
  // way, sorting the encoded endpoint sorts in order of x1, y1, x2, y2
  // efficiently.
  std::vector<std::pair<uint64_t, int> > sorted_endpoints;
  sorted_endpoints.reserve(features1_.keypoints.size());
  for (int i = 0; i < features1_.keypoints.size(); i++) {
    if (ContainsKey(matched_features1_, i)) {
      continue;
    }

    const Eigen::Vector2d point(features1_.keypoints[i].x(),
                                features1_.keypoints[i].y());
    // Compute the epipolar line.
    Eigen::Vector3d epipolar_line = fundamental_matrix * point.homogeneous();
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

    sorted_endpoints.emplace_back(EncodeLineEndpoints(line_endpoints), i);
  }

  // Sort the endpoints.
  std::sort(sorted_endpoints.begin(), sorted_endpoints.end());

  // Cluster the epilines basd on the epipolar line endpoints.
  const double sq_max_distance_pixels =
      options_.guided_matching_max_distance_pixels *
      options_.guided_matching_max_distance_pixels;

  for (int i = 0; i < sorted_endpoints.size(); i++) {
    std::vector<Eigen::Vector2d> line_endpoints;
    DecodeLineEndpoints(sorted_endpoints[i].first, &line_endpoints);

    // See if an epipolar endpoint exists nearby. If it does not, then add a new
    // endpoint group.
    if (i == 0 ||
        (epiline_groups->back().endpoints[0] - line_endpoints[0])
                .squaredNorm() > sq_max_distance_pixels) {
      EpilineGroup epiline_group;
      epiline_group.endpoints = line_endpoints;
      epiline_groups->emplace_back(epiline_group);
    }

    // Assign the feature to the most recent epipolar line.
    epiline_groups->back().features.emplace_back(sorted_endpoints[i].second);
  }
  LOG(INFO) << "Created " << epiline_groups->size() << " groups from "
            << sorted_endpoints.size() << " original epilines.";
}

void GuidedEpipolarMatcher::FindFeaturesNearEpipolarLines(
    const EpilineGroup& epiline_group,
    std::vector<int>* candidate_keypoint_indices) {
  static const int kMinNumMatchesFound = 25;

  const std::vector<Eigen::Vector2d>& line_endpoints = epiline_group.endpoints;

  // The number of steps required to "walk" between the line endpoints.
  const int num_steps =
    static_cast<int>((line_endpoints[1] - line_endpoints[0]).norm() /
                     options_.guided_matching_max_distance_pixels);

  // Sample the epipolar line equally between the points where it intersects
  // the features bounding box.
  std::unordered_set<int> candidate_keypoints;
  const Eigen::Vector2d line_delta =
      (line_endpoints[0] - line_endpoints[1]) / static_cast<double>(num_steps);
  Eigen::Vector2d sample_point = line_endpoints[1];

  for (int i = 0; i < num_steps; i++) {
    sample_point += line_delta;

    // Find the cell center among all grids that is closest and add the
    // keypoints belonging to that cell.
    std::vector<int> new_keypoints;
    FindClosestCellAndKeypoints(sample_point, &new_keypoints);
    candidate_keypoints.insert(new_keypoints.begin(), new_keypoints.end());
  }

  // If we do not have enough features then the lowes ratio test is not
  // meaningful. Add some random features here so that lowes ratio is more
  // informative of whether we have a good match or not.
  if (candidate_keypoints.size() < kMinNumMatchesFound) {
    const int num_current_keypoints = candidate_keypoints.size();
    for (int i = num_current_keypoints; i < kMinNumMatchesFound; i++) {
      candidate_keypoints.insert(RandInt(0, features2_.keypoints.size() - 1));
    }
  }
  candidate_keypoint_indices->insert(candidate_keypoint_indices->end(),
                                     candidate_keypoints.begin(),
                                     candidate_keypoints.end());
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
    const std::vector<int>& query_feature_indices,
    const std::vector<int>& candidate_feature_indices,
    std::vector<std::vector<float> >* nn_distances,
    std::vector<std::vector<int> >* nn_indices) {
  static const int kNumNearestNeighbors = 2;
  static const int kMinNumLeafsVisited = 50;
  const int num_descriptor_dimensions = features1_.descriptors[0].size();

  // Gather the query descriptors.
  Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
      query_descriptors(query_feature_indices.size(),
                        num_descriptor_dimensions);
  for (int i = 0; i < query_feature_indices.size(); i++) {
    const int match_index = query_feature_indices[i];
    query_descriptors.row(i) = features1_.descriptors[match_index];
  }
  flann::Matrix<float> flann_query_descriptors(query_descriptors.data(),
                                               query_descriptors.rows(),
                                               query_descriptors.cols());

  // Gather the candidate matching descriptors.
  Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
      candidate_descriptors(candidate_feature_indices.size(),
                            num_descriptor_dimensions);
  for (int i = 0; i < candidate_feature_indices.size(); i++) {
    const int match_index = candidate_feature_indices[i];
    candidate_descriptors.row(i) = features2_.descriptors[match_index];
  }

  // Create the searchable KD-tree with FLANN.
  flann::Matrix<float> flann_candidate_descriptors(
      candidate_descriptors.data(), candidate_descriptors.rows(),
      candidate_descriptors.cols());

  flann::Index<flann::L2<float> > flann_kd_tree(
      flann_candidate_descriptors, flann::KDTreeSingleIndexParams());
  flann_kd_tree.buildIndex();

  // Query the KD-tree to get the top 2 nearest neighbors.
  const int max_leafs_to_check =
      std::max(static_cast<int>(candidate_descriptors.rows() * 0.2),
               kMinNumLeafsVisited);
  flann_kd_tree.knnSearch(flann_query_descriptors, *nn_indices, *nn_distances,
                          kNumNearestNeighbors,
                          flann::SearchParams(max_leafs_to_check));

  // Output the top 2 matches.
  for (int i = 0; i < query_feature_indices.size(); i++) {
    // Change the NN indices to be the feature index instead of the index of the
    // FLANN matrix.
    for (int j = 0; j < kNumNearestNeighbors; j++) {
      const int flann_index = (*nn_indices)[i][j];
      (*nn_indices)[i][j] = candidate_feature_indices[flann_index];
    }
  }
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
