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

#include "theia/sfm/pose/estimate_positions_linear.h"

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseQR>
#include <glog/logging.h>
#include <unordered_map>
#include <vector>

#include "theia/util/map_util.h"
#include "theia/util/threadpool.h"
#include "theia/sfm/reconstruction.h"
#include "theia/sfm/types.h"
#include "theia/sfm/view_triplet.h"

namespace theia {

using Eigen::Matrix3d;
using Eigen::Vector3d;

namespace {

// Adds the 3x3 matrix to the triplet list with the top-left corner in the
// sparse matrix corresponding to (row, col).
void Add3x3MatrixToTripletList(
    const Matrix3d& mat,
    const int row,
    const int col,
    std::vector<Eigen::Triplet<double> >* triplet_list) {
  triplet_list->emplace_back(row + 0, col + 0, mat(0, 0));
  triplet_list->emplace_back(row + 1, col + 0, mat(1, 0));
  triplet_list->emplace_back(row + 2, col + 0, mat(2, 0));
  triplet_list->emplace_back(row + 0, col + 1, mat(0, 1));
  triplet_list->emplace_back(row + 1, col + 1, mat(1, 1));
  triplet_list->emplace_back(row + 2, col + 1, mat(2, 1));
  triplet_list->emplace_back(row + 0, col + 2, mat(0, 2));
  triplet_list->emplace_back(row + 1, col + 2, mat(1, 2));
  triplet_list->emplace_back(row + 2, col + 2, mat(2, 2));
}

inline Matrix3d AngleAxisToRotationMatrix(const Vector3d angle_axis) {
  const double angle = angle_axis.norm();
  const Eigen::AngleAxisd rotation_aa(angle, angle_axis / angle);
  return rotation_aa.toRotationMatrix();
}

// Adds a triplet constraint to the linear system. The weight of the constraint
// (w), the global orientations, baseline (ratios), and view triplet information
// are needed to form the constraint.
void AddTripletConstraint(const ViewTriplet& triplet,
                          const int start_row,
                          const std::vector<int> cols,
                          const double w,
                          const std::vector<Vector3d>& orientations,
                          const Vector3d& baselines,
                          std::vector<Eigen::Triplet<double> >* triplet_list) {
  // Relative camera positions.
  const Matrix3d orientation0 = AngleAxisToRotationMatrix(orientations[0]);
  const Matrix3d orientation1 = AngleAxisToRotationMatrix(orientations[1]);
  const Vector3d& t01 =
      -orientation0.transpose() * triplet.info_one_two.position_2;
  const Vector3d& t02 =
      -orientation0.transpose() * triplet.info_one_three.position_2;
  const Vector3d& t12 =
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
  Matrix3d m1 =
      (-s_201 * r201 + r012.transpose() / s_012 + Matrix3d::Identity()) * w;
  Matrix3d m2 =
      (s_201 * r201 - r012.transpose() / s_012 + Matrix3d::Identity()) * w;
  Matrix3d m3 = -2.0 * w * Matrix3d::Identity();
  Add3x3MatrixToTripletList(m1, start_row, cols[0], triplet_list);
  Add3x3MatrixToTripletList(m2, start_row, cols[1], triplet_list);
  Add3x3MatrixToTripletList(m3, start_row, cols[2], triplet_list);

  // Assume t02 is perfect and solve for c1.
  m1 = (-r201.transpose() / s_201 + s_120 * r120 + Matrix3d::Identity()) * w;
  m2 = -2.0 * w * Matrix3d::Identity();
  m3 = (r201.transpose() / s_201 - s_120 * r120 + Matrix3d::Identity()) * w;
  Add3x3MatrixToTripletList(m1, start_row + 3, cols[0], triplet_list);
  Add3x3MatrixToTripletList(m2, start_row + 3, cols[1], triplet_list);
  Add3x3MatrixToTripletList(m3, start_row + 3, cols[2], triplet_list);

  // Assume t12 is perfect and solve for c0.
  m1 = -2.0  * w * Matrix3d::Identity();
  m2 = (-s_012 * r012 + r120.transpose() / s_120 + Matrix3d::Identity()) * w;
  m3 = (s_012 * r012 - r120.transpose() / s_120 + Matrix3d::Identity()) * w;
  Add3x3MatrixToTripletList(m1, start_row + 6, cols[0], triplet_list);
  Add3x3MatrixToTripletList(m2, start_row + 6, cols[1], triplet_list);
  Add3x3MatrixToTripletList(m3, start_row + 6, cols[2], triplet_list);
}

}  // namespace

LinearPositionEstimator::LinearPositionEstimator(
    const LinearPositionEstimatorOptions& options,
    const Reconstruction& reconstruction)
    : options_(options), reconstruction_(reconstruction) {
  CHECK_GT(options.num_threads, 0);
  CHECK_GT(options.max_translation_angle_change_degrees, 0);
  CHECK_GT(options.max_reprojection_error_pixels, 0);
  CHECK_GT(options.min_good_reprojections, 0);
}

bool LinearPositionEstimator::EstimatePositions(
    const std::vector<ViewTriplet>& triplets,
    const std::unordered_map<ViewId, Vector3d>& orientations,
    std::unordered_map<ViewId, Vector3d>* positions) {
  CHECK_NOTNULL(positions)->clear();

  // Count the number of times each view is in a triplet.
  for (int i = 0; i < triplets.size(); i++) {
    num_triplets_for_view_[triplets[i].view_ids[0]] =
        num_triplets_for_view_[triplets[i].view_ids[0]] + 1;
    num_triplets_for_view_[triplets[i].view_ids[1]] =
        num_triplets_for_view_[triplets[i].view_ids[1]] + 1;
    num_triplets_for_view_[triplets[i].view_ids[2]] =
        num_triplets_for_view_[triplets[i].view_ids[2]] + 1;

    // Determine the order of the views in the linear system.
    InsertIfNotPresent(&linear_system_index_,
                       triplets[i].view_ids[0],
                       linear_system_index_.size());
    InsertIfNotPresent(&linear_system_index_,
                       triplets[i].view_ids[1],
                       linear_system_index_.size());
    InsertIfNotPresent(&linear_system_index_,
                       triplets[i].view_ids[2],
                       linear_system_index_.size());
  }

  // Baselines where (x, y, z) corresponds to the baseline of the first, second,
  // and third view pair in the triplet.
  std::vector<Vector3d> baselines;
  ComputeBaselineRatios(triplets, &baselines);

  Eigen::SparseMatrix<double> constraint_matrix;
  CreateLinearSystem(triplets, orientations, baselines, &constraint_matrix);

  const Eigen::SparseMatrix<double> aTa =
      constraint_matrix.transpose() * constraint_matrix;

  // TODO(cmsweeney): Solve the eigenvalue problem.
  const Eigen::VectorXd& solution;

  // Add the solutions to the output.
  for (int i = 0; i < solution.size() / 3; i++) {
    const ViewId view_id = FindOrDie(linear_system_index_, i);
    (*positions)[view_id] = solution.segment<3>(i * 3);
  }

  return true;
}

void LinearPositionEstimator::ComputeBaselineRatios(
    const std::vector<ViewTriplet>& triplets,
    std::vector<Vector3d>* baselines) {
  CHECK_NOTNULL(baselines)->resize(triplets.size());

  ThreadPool pool(options_.num_threads);
  for (int i = 0; i < triplets.size(); i++) {
    pool.Add(&ComputeBaselineRatioForTriplet,
             this,
             triplets[i],
             &baselines->at(i));
  }
}

void LinearPositionEstimator::ComputeBaselineForTriplet(
    const ViewTriplet& triplet, Vector3d* baseline) {
  baseline->setZero();

  const View& view1 = *reconstruction.View(triplet.view_ids[0]);
  const View& view2 = *reconstruction.View(triplet.view_ids[1]);
  const View& view3 = *reconstruction.View(triplet.view_ids[2]);

  // Find common tracks.
  std::vector<TrackId> common_tracks =
      FindCommonTracksInTriplet(triplet.view_ids);

  // Compute the baseline ratios for all tracks.
  Eigen::Vector4d point12, point13, point23;
  double depth1_12, depth2_12, depth1_13, depth3_13, depth2_23, depth3_23;
  for (const TrackId track_id : common_tracks) {
    const Feature& feature1 = GetNormalizedFeature(view1, track_id);
    const Feature& feature2 = GetNormalizedFeature(view2, track_id);
    const Feature& feature3 = GetNormalizedFeature(view3, track_id);

    // Compute triangulation from views 1, 2.
    GetTriangulatedPointDepths(info_one_two,
                               feature1, feature2,
                               &depth1_12, &depth2_12);

    // Compute triangulation from views 1, 3.
    GetTriangulatedPointDepths(info_one_three,
                               feature1, feature3,
                               &depth1_13, &depth3_13);

    // Compute triangulation from views 2, 3.
    GetTriangulatedPointDepths(info_two_three,
                               feature2, feature3,
                               &depth2_23, &depth3_23);
    *baseline +=
        Vector3d(1.0, depth1_13 / depth1_12, depth2_23 / depth2_12);
  }

  // Take the mean as the baseline ratios.
  *baseline /= static_cast<double>(common_tracks.size());
}

std::vector<TrackId> LinearPositionEstimator::FindCommonTracks(
    const ViewId view_ids[3]) {
  // Get all track ids for the views.
  std::vector<TrackId> view1_tracks =
      reconstruction_.View(view_ids[0])->TrackIds();
  std::vector<TrackId> view2_tracks =
      reconstruction_.View(view_ids[1])->TrackIds();
  std::vector<TrackId> view3_tracks =
      reconstruction_.View(view_ids[2])->TrackIds();

  // Sort the track id.
  std::sort(view1_tracks.begin(), view1_tracks.end());
  std::sort(view2_tracks.begin(), view2_tracks.end());
  std::sort(view3_tracks.begin(), view3_tracks.end());

  // Get the intersection of all the tracks.
  std::vector<TrackId> view12_intersection;
  std::intersection(view1_tracks.begin(), view1_tracks.end(),
                    view2_tracks.begin(), view2_tracks.end(),
                    std::back_inserter(view12_intersection));
  std::vector<TrackId> view123_intersection;
  std::intersection(view12_intersection.begin(), view12_intersection.end(),
                    view3_tracks.begin(), view3_tracks.end(),
                    std::back_inserter(view123_intersection));
  return view123_intersection;
}

void LinearPositionEstimator::GetTriangulatedPointDepths(
    const TwoViewInfo& info,
    const Feature& feature1,
    const Feature& feature2,
    double* depth1,
    double* depth2) {
  // TODO: Normalize the points.
  const Feature normalized_feature1, normalized_feature2;

  // Triangulate point.
  Eigen::Vector4d point;
  const Vector3d origin1, origin2, ray_direction1, ray_direction2;
  Triangulate(origin1, ray_direction1, origin2, ray_direction2, &point);

  // Compute depths.
  const Vector3d point3d = point.hnormalized();
  *depth1 = point3d.norm();
  *depth2 = (point3d - info.position_2).norm();
}

// Sets up the linear system with the constraints that each triplet adds.
void LinearPositionEstimator::CreateLinearSystem(
    const std::vector<ViewTriplet>& triplets,
    const std::unordered_map<ViewId, Vector3d>& orientations,
    std::vector<Vector3d>& baselines,
    Eigen::SparseMatrix<double>* constraint_matrix) {
  const int num_views = num_triplets_for_view_.size();

  constraint_matrix->resize(triplets.size() * 9, num_views * 3);

  std::vector<Eigen::Triplet<double> > triplet_list;

  for (int i = 0; i < triplets.size(); i++) {
    const ViewTriplet& triplet = triplets[i];
    const std::vector<Vector3d> triplet_orientations = {
      FindOrDie(orientations, triplet.view_ids[0]),
      FindOrDie(orientations, triplet.view_ids[1]),
      FindOrDie(orientations, triplet.view_ids[2])
    };
    // Get the row and columns that we will modify.
    const std::vector<int> cols = {
        static_cast<int>(3 *
                         FindOrDie(linear_system_index_, triplet.view_ids[0])),
        static_cast<int>(3 *
                         FindOrDie(linear_system_index_, triplet.view_ids[1])),
        static_cast<int>(3 *
                         FindOrDie(linear_system_index_, triplet.view_ids[2]))};

    const double w =
        1.0 / sqrt(std::min({num_triplets_for_view_[triplet.view_ids[0]],
                             num_triplets_for_view_[triplet.view_ids[1]],
                             num_triplets_for_view_[triplet.view_ids[2]]}));

    AddTripletConstraint(triplet,
                         9 * i,
                         cols,
                         w,
                         triplet_orientations,
                         baselines[i],
                         &triplet_list);
  }
  constraint_matrix->setFromTriplets(triplet_list.begin(), triplet_list.end());
}

}  // namespace theia
