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

#include "theia/sfm/global_pose_estimation/linear_position_estimator.h"

#include <ceres/rotation.h>
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <glog/logging.h>
#include <algorithm>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "spectra/include/SymEigsSolver.h"

#include "theia/math/matrix/linear_operator.h"
#include "theia/sfm/find_common_tracks_in_views.h"
#include "theia/sfm/global_pose_estimation/compute_triplet_baseline_ratios.h"
#include "theia/sfm/reconstruction.h"
#include "theia/sfm/types.h"
#include "theia/sfm/view_graph/triplet_extractor.h"
#include "theia/sfm/view_triplet.h"
#include "theia/util/map_util.h"
#include "theia/util/threadpool.h"

namespace theia {

using Eigen::Matrix3d;
using Eigen::Vector3d;

namespace {

// Adds the constraint from the triplet to the symmetric matrix. Our standard
// constraint matrix A is a 3M x 3N matrix with M triplet constraints and N
// cameras. We seek to construct A^t * A directly. For each triplet constraint
// in our matrix A (i.e. a 3-row block), we can compute the corresponding
// entries in A^t * A with the following summation:
//
//   A^t * A += Row(i)^t * Row(i)
//
// for each triplet constraint i.
void AddTripletConstraintToSymmetricMatrix(
    const std::vector<Matrix3d>& constraints,
    const std::vector<int>& view_indices,
    std::unordered_map<std::pair<int, int>, double>* sparse_matrix_entries) {
  // Construct Row(i)^t * Row(i). If we denote the row as a block matrix:
  //
  //   Row(i) = [A | B | C]
  //
  // then we have:
  //
  //   Row(i)^t * Row(i) = [A | B | C]^t * [A | B | C]
  //                     = [ A^t * A  |  A^t * B  |  A^t * C]
  //                       [ B^t * A  |  B^t * B  |  B^t * C]
  //                       [ C^t * A  |  C^t * B  |  C^t * C]
  //
  // Since A^t * A is symmetric, we only store the upper triangular portion.
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      // Skip any block entries that correspond to the lower triangular portion
      // of the matrix.
      if (view_indices[i] > view_indices[j]) {
        continue;
      }

      // Compute the A^t * B, etc. matrix.
      const Eigen::Matrix3d symmetric_constraint =
          constraints[i].transpose() * constraints[j];

      // Add to the 3x3 block corresponding to (i, j)
      for (int r = 0; r < 3; r++) {
        for (int c = 0; c < 3; c++) {
          const std::pair<int, int> row_col(view_indices[i] + r,
                                            view_indices[j] + c);
          (*sparse_matrix_entries)[row_col] += symmetric_constraint(r, c);
        }
      }
    }
  }
}

inline Matrix3d AngleAxisToRotationMatrix(const Vector3d angle_axis) {
  const double angle = angle_axis.norm();
  const Eigen::AngleAxisd rotation_aa(angle, angle_axis / angle);
  return rotation_aa.toRotationMatrix();
}

// Adds a triplet constraint to the linear system. The weight of the constraint
// (w), the global orientations, baseline (ratios), and view triplet information
// are needed to form the constraint.
void AddTripletConstraintToSparseMatrix(
    const ViewTriplet& triplet,
    const std::vector<int>& view_indices,
    const double w,
    const std::vector<Vector3d>& orientations,
    const Vector3d& baselines,
    std::unordered_map<std::pair<int, int>, double>* sparse_matrix_entries) {
  // Relative camera positions.
  const Matrix3d orientation0 = AngleAxisToRotationMatrix(orientations[0]);
  const Matrix3d orientation1 = AngleAxisToRotationMatrix(orientations[1]);
  const Vector3d t01 =
      -orientation0.transpose() * triplet.info_one_two.position_2;
  const Vector3d t02 =
      -orientation0.transpose() * triplet.info_one_three.position_2;
  const Vector3d t12 =
      -orientation1.transpose() * triplet.info_two_three.position_2;

  // Rotations between the translation vectors.
  const Matrix3d r012 =
      Eigen::Quaterniond::FromTwoVectors(t12, -t01).toRotationMatrix();
  const Matrix3d r201 =
      Eigen::Quaterniond::FromTwoVectors(t01, t02).toRotationMatrix();
  const Matrix3d r120 =
      Eigen::Quaterniond::FromTwoVectors(-t02, -t12).toRotationMatrix();

  // Baselines ratios.
  const double s_012 = baselines[0] / baselines[2];
  const double s_201 = baselines[1] / baselines[0];
  const double s_120 = baselines[2] / baselines[1];

  // Assume that t01 is perfect and solve for c2.
  std::vector<Eigen::Matrix3d> constraints(3);
  constraints[0] =
      (-s_201 * r201 + r012.transpose() / s_012 + Matrix3d::Identity()) * w;
  constraints[1] =
      (s_201 * r201 - r012.transpose() / s_012 + Matrix3d::Identity()) * w;
  constraints[2] = -2.0 * w * Matrix3d::Identity();
  AddTripletConstraintToSymmetricMatrix(constraints,
                                        view_indices,
                                        sparse_matrix_entries);

  // Assume t02 is perfect and solve for c1.
  constraints[0] =
      (-r201.transpose() / s_201 + s_120 * r120 + Matrix3d::Identity()) * w;
  constraints[1] = -2.0 * w * Matrix3d::Identity();
  constraints[2] =
      (r201.transpose() / s_201 - s_120 * r120 + Matrix3d::Identity()) * w;
  AddTripletConstraintToSymmetricMatrix(constraints,
                                        view_indices,
                                        sparse_matrix_entries);

  // Assume t12 is perfect and solve for c0.
  constraints[0] = -2.0  * w * Matrix3d::Identity();
  constraints[1] =
      (-s_012 * r012 + r120.transpose() / s_120 + Matrix3d::Identity()) * w;
  constraints[2] =
      (s_012 * r012 - r120.transpose() / s_120 + Matrix3d::Identity()) * w;
  AddTripletConstraintToSymmetricMatrix(constraints,
                                        view_indices,
                                        sparse_matrix_entries);
}

// Returns true if the vector R1 * (c2 - c1) is in the same direction as t_12.
bool VectorsAreSameDirection(const Vector3d& position1,
                             const Vector3d& position2,
                             const Vector3d& rotation1,
                             const Vector3d& relative_position12) {
  const Vector3d global_relative_position =
      (position2 - position1).normalized();
  Vector3d rotated_relative_position;
  ceres::AngleAxisRotatePoint(rotation1.data(),
                              global_relative_position.data(),
                              rotated_relative_position.data());
  return rotated_relative_position.dot(relative_position12) > 0;
}

}  // namespace

LinearPositionEstimator::LinearPositionEstimator(
    const Options& options,
    const Reconstruction& reconstruction)
    : options_(options), reconstruction_(reconstruction) {
  CHECK_GT(options.num_threads, 0);
}

bool LinearPositionEstimator::EstimatePositions(
    const std::unordered_map<ViewIdPair, TwoViewInfo>& view_pairs,
    const std::unordered_map<ViewId, Vector3d>& orientations,
    std::unordered_map<ViewId, Vector3d>* positions) {
  CHECK_NOTNULL(positions)->clear();

  // Extract triplets from the view pairs. As of now, we only consider the
  // largest connected triplet in the viewing graph.
  // TODO(cmsweeney): Utilize all connected triplet graphs.
  VLOG(2) << "Extracting triplets from the viewing graph.";
  TripletExtractor extractor;
  std::vector<std::vector<ViewTriplet> > triplets_vec;
  CHECK(extractor.ExtractTripletsFromViewPairs(view_pairs, &triplets_vec));
  // Find the largest triplet.
  int largest_triplet_graph = 0;
  for (int i = 1; i < triplets_vec.size(); i++) {
    if (triplets_vec[i].size() > triplets_vec[largest_triplet_graph].size()) {
      largest_triplet_graph = i;
    }
  }

  // Add all triplets in the largest connected component to the problem.
  for (int i = 0; i < triplets_vec[largest_triplet_graph].size(); i++) {
    AddTripletConstraint(triplets_vec[largest_triplet_graph][i]);
  }

  return EstimatePositions(orientations, positions);
}

// An alternative interface is to instead add triplets one by one to linear
// estimator. This allows for adding redundant observations of triplets, which
// may be useful if there are multiple estimates of the data.
void LinearPositionEstimator::AddTripletConstraint(
    const ViewTriplet& view_triplet) {
  triplets_.emplace_back(view_triplet);
  num_triplets_for_view_[view_triplet.view_ids[0]] += 1;
  num_triplets_for_view_[view_triplet.view_ids[1]] += 1;
  num_triplets_for_view_[view_triplet.view_ids[2]] += 1;

  // Determine the order of the views in the linear system. We subtract 1 from
  // the linear system index so that the first position added to the system
  // will be set constant (index of -1 is intentionally not evaluated later).
  InsertIfNotPresent(&linear_system_index_,
                     view_triplet.view_ids[0],
                     linear_system_index_.size() - 1);
  InsertIfNotPresent(&linear_system_index_,
                     view_triplet.view_ids[1],
                     linear_system_index_.size() - 1);
  InsertIfNotPresent(&linear_system_index_,
                     view_triplet.view_ids[2],
                     linear_system_index_.size() - 1);
}

bool LinearPositionEstimator::EstimatePositions(
    const std::unordered_map<ViewId, Eigen::Vector3d>& orientations,
    std::unordered_map<ViewId, Eigen::Vector3d>* positions) {
  VLOG(2) << "Determining baseline ratios within each triplet...";
  ComputeBaselineRatios();

  VLOG(2) << "Building the constraint matrix...";
  // Create the linear system based on triplet constraints.
  Eigen::SparseMatrix<double> constraint_matrix;
  CreateLinearSystem(orientations, &constraint_matrix);

  // Solve for positions by examining the smallest eigenvalues. Since we have
  // set one position constant at the origin, we only need to solve for the
  // eigenvector corresponding to the smallest eigenvalue. This can be done
  // efficiently with inverse power iterations.
  VLOG(2) << "Solving for positions from the sparse eigenvalue problem...";
  SparseSymShiftSolveLLT op(constraint_matrix);
  Spectra::SymEigsShiftSolver<double, Spectra::LARGEST_MAGN,
                              SparseSymShiftSolveLLT> eigs(&op, 1, 6, 0.0);
  eigs.init();
  eigs.compute();

  // Compute with power iterations.
  const Eigen::VectorXd solution = eigs.eigenvectors().col(0);

  // Add the solutions to the output. Set the position with an index of -1 to
  // be at the origin.
  for (const auto& view_index : linear_system_index_) {
    if (view_index.second < 0) {
      (*positions)[view_index.first].setZero();
    } else {
      (*positions)[view_index.first] =
          solution.segment<3>(view_index.second * 3);
    }
  }

  // Flip the sign of the positions if necessary.
  FlipSignOfPositionsIfNecessary(orientations, positions);

  return true;
}

void LinearPositionEstimator::ComputeBaselineRatios() {
  baselines_.resize(triplets_.size());

  // Only estimate the baseline ratio from features if requested. Otherwise, use
  // the scale of the relative translations provided from the input as the
  // baselines.
  if (!options_.estimate_relative_baselines_from_features) {
    for (int i = 0; i < triplets_.size(); i++) {
      // NOTE: We do not need to normalize the relative positions here since
      // they are only used to compute the rotation between relative
      // translations in the triplet constraint (i.e. scale is not considered).
      baselines_[i][0] = triplets_[i].info_one_two.position_2.norm();
      baselines_[i][1] = triplets_[i].info_one_three.position_2.norm();
      baselines_[i][2] = triplets_[i].info_two_three.position_2.norm();
    }
  } else {
    // Baselines where (x, y, z) corresponds to the baseline of the first,
    // second, and third view pair in the triplet.
    ThreadPool pool(options_.num_threads);
    for (int i = 0; i < triplets_.size(); i++) {
      pool.Add(&LinearPositionEstimator::ComputeBaselineRatioForTriplet,
               this,
               triplets_[i],
               &baselines_[i]);
    }
  }
}

void LinearPositionEstimator::ComputeBaselineRatioForTriplet(
    const ViewTriplet& triplet, Vector3d* baseline) {
  baseline->setZero();

  const View& view1 = *reconstruction_.View(triplet.view_ids[0]);
  const View& view2 = *reconstruction_.View(triplet.view_ids[1]);
  const View& view3 = *reconstruction_.View(triplet.view_ids[2]);

  // Find common tracks.
  const std::vector<ViewId> triplet_view_ids = {
      triplet.view_ids[0], triplet.view_ids[1], triplet.view_ids[2]};
  const std::vector<TrackId>& common_tracks =
      FindCommonTracksInViews(reconstruction_, triplet_view_ids);

  // Normalize all features.
  std::vector<Feature> feature1, feature2, feature3;
  feature1.reserve(common_tracks.size());
  feature2.reserve(common_tracks.size());
  feature3.reserve(common_tracks.size());
  for (const TrackId track_id : common_tracks) {
    feature1.emplace_back(GetNormalizedFeature(view1, track_id));
    feature2.emplace_back(GetNormalizedFeature(view2, track_id));
    feature3.emplace_back(GetNormalizedFeature(view3, track_id));
  }

  // Get the baseline ratios.
  ComputeTripletBaselineRatios(triplet, feature1, feature2, feature3, baseline);
}

// Sets up the linear system with the constraints that each triplet adds.
void LinearPositionEstimator::CreateLinearSystem(
    const std::unordered_map<ViewId, Vector3d>& orientations,
    Eigen::SparseMatrix<double>* constraint_matrix) {
  const int num_views = num_triplets_for_view_.size();

  std::unordered_map<std::pair<int, int>, double> sparse_matrix_entries;
  sparse_matrix_entries.reserve(27 * num_triplets_for_view_.size());
  for (int i = 0; i < triplets_.size(); i++) {
    const ViewTriplet& triplet = triplets_[i];
    const std::vector<Vector3d> triplet_orientations = {
      FindOrDie(orientations, triplet.view_ids[0]),
      FindOrDie(orientations, triplet.view_ids[1]),
      FindOrDie(orientations, triplet.view_ids[2])
    };
    // Get the index of each camera in the sparse matrix.
    const std::vector<int> view_indices = {
        static_cast<int>(3 *
                         FindOrDie(linear_system_index_, triplet.view_ids[0])),
        static_cast<int>(3 *
                         FindOrDie(linear_system_index_, triplet.view_ids[1])),
        static_cast<int>(3 *
                         FindOrDie(linear_system_index_, triplet.view_ids[2]))};

    const double w =
        1.0 /
        std::sqrt(std::min({num_triplets_for_view_[triplet.view_ids[0]],
                            num_triplets_for_view_[triplet.view_ids[1]],
                            num_triplets_for_view_[triplet.view_ids[2]]}));
    AddTripletConstraintToSparseMatrix(triplet,
                                       view_indices,
                                       w,
                                       triplet_orientations,
                                       baselines_[i],
                                       &sparse_matrix_entries);
  }

  // Set the sparse matrix from the container of the accumulated entries.
  std::vector<Eigen::Triplet<double> > triplet_list;
  triplet_list.reserve(sparse_matrix_entries.size());
  for (const auto& sparse_matrix_entry : sparse_matrix_entries) {
    // Skip this entry if the indices are invalid. This only occurs when we
    // encounter a constraint with the constant camera (which has a view index
    // of -1).
    if (sparse_matrix_entry.first.first < 0 ||
        sparse_matrix_entry.first.second < 0) {
      continue;
    }
    triplet_list.emplace_back(sparse_matrix_entry.first.first,
                              sparse_matrix_entry.first.second,
                              sparse_matrix_entry.second);
  }

  // We construct the constraint matrix A^t * A directly, which is an
  // N - 1 x N - 1 matrix where N is the number of cameras (and 3 entries per
  // camera, corresponding to the camera position entries).
  constraint_matrix->resize((num_views - 1) * 3, (num_views - 1) * 3);
  constraint_matrix->setFromTriplets(triplet_list.begin(), triplet_list.end());
}

Feature LinearPositionEstimator::GetNormalizedFeature(const View& view,
                                                      const TrackId track_id) {
  Feature feature = *view.GetFeature(track_id);
  const Camera& camera = view.Camera();
  Eigen::Vector3d normalized_feature =
      camera.PixelToNormalizedCoordinates(feature);
  return normalized_feature.hnormalized();
}

void LinearPositionEstimator::FlipSignOfPositionsIfNecessary(
    const std::unordered_map<ViewId, Vector3d>& orientation,
    std::unordered_map<ViewId, Vector3d>* positions) {
  // If this value is below zero, then we should flip the sign.
  int correct_sign_votes = 0;

  std::unordered_set<ViewIdPair> pairs_visited;
  for (const ViewTriplet& triplet : triplets_) {
    const ViewIdPair id_pair_12(triplet.view_ids[0], triplet.view_ids[1]);
    const ViewIdPair id_pair_13(triplet.view_ids[0], triplet.view_ids[2]);
    const ViewIdPair id_pair_23(triplet.view_ids[1], triplet.view_ids[2]);

    // It is not guaranteed that these positions are needed, but it should be
    // very fast to fetch them regardless so this will not be a significant
    // performance cost.
    const Vector3d& position1 = FindOrDie(*positions, triplet.view_ids[0]);
    const Vector3d& position2 = FindOrDie(*positions, triplet.view_ids[1]);
    const Vector3d& position3 = FindOrDie(*positions, triplet.view_ids[2]);

    // Check the relative translation of views 1 and 2 in the triplet.
    if (!ContainsKey(pairs_visited, id_pair_12)) {
      pairs_visited.insert(id_pair_12);
      correct_sign_votes +=
          (VectorsAreSameDirection(position1, position2,
                                   FindOrDie(orientation, id_pair_12.first),
                                   triplet.info_one_two.position_2))
              ? 1
              : -1;
    }

    // Check the relative translation of views 1 and 3 in the triplet.
    if (!ContainsKey(pairs_visited, id_pair_13)) {
      pairs_visited.insert(id_pair_13);
      correct_sign_votes +=
          (VectorsAreSameDirection(position1, position3,
                                   FindOrDie(orientation, id_pair_13.first),
                                   triplet.info_one_three.position_2))
              ? 1
              : -1;
    }

    // Check the relative translation of views 2 and 3 in the triplet.
    if (!ContainsKey(pairs_visited, id_pair_23)) {
      pairs_visited.insert(id_pair_23);
      correct_sign_votes +=
          (VectorsAreSameDirection(position2, position3,
                                   FindOrDie(orientation, id_pair_23.first),
                                   triplet.info_two_three.position_2))
              ? 1
              : -1;
    }
  }

  // If the sign of the votes is below zero, we must flip the sign of all
  // position estimates.
  if (correct_sign_votes < 0) {
    const int num_correct_votes =
        (pairs_visited.size() + correct_sign_votes) / 2;
    VLOG(2) << "Sign of the positions was incorrect: " << num_correct_votes
            << " of " << pairs_visited.size()
            << " relative translations had the correct sign. "
               "Flipping the sign of the camera positions.";
    for (auto& position : *positions) {
      position.second *= -1.0;
    }
  }
}

}  // namespace theia
