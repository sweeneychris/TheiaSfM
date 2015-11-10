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

#include "theia/sfm/extract_maximally_parallel_rigid_subgraph.h"

#include <ceres/rotation.h>
#include <Eigen/Core>
#include <Eigen/LU>

#include <algorithm>
#include <unordered_map>
#include <vector>

#include "theia/sfm/pose/util.h"
#include "theia/sfm/types.h"
#include "theia/sfm/view_graph/view_graph.h"
#include "theia/util/map_util.h"

namespace theia {
namespace {

void FormAngleMeasurementMatrix(
    const std::unordered_map<ViewId, Eigen::Vector3d>& orientations,
    const ViewGraph& view_graph,
    const std::unordered_map<ViewId, int>& view_ids_to_index,
    Eigen::MatrixXd* angle_measurements) {
  const auto& view_pairs = view_graph.GetAllEdges();
  angle_measurements->setZero();

  // Set up the matrix such that t_{i,j} x (c_j - c_i) = 0.
  int i = 0;
  for (const auto& view_pair : view_pairs) {
    // Get t_{i,j} and rotate it such that it is oriented in the global
    // reference frame.
    Eigen::Matrix3d world_to_view1_rotation;
    ceres::AngleAxisToRotationMatrix(
        FindOrDie(orientations, view_pair.first.first).data(),
        ceres::ColumnMajorAdapter3x3(world_to_view1_rotation.data()));
    const Eigen::Vector3d rotated_translation =
        world_to_view1_rotation.transpose() * view_pair.second.position_2;
    const Eigen::Matrix3d cross_product_mat =
        CrossProductMatrix(rotated_translation);

    // Find the column locations of the two views.
    const int view1_col =
        3 * FindOrDie(view_ids_to_index, view_pair.first.first);
    const int view2_col =
        3 * FindOrDie(view_ids_to_index, view_pair.first.second);

    angle_measurements->block<3, 3>(3 * i, view1_col) = -cross_product_mat;
    angle_measurements->block<3, 3>(3 * i, view2_col) = cross_product_mat;
    ++i;
  }
}

// Computes the cosine distance in each dimension x, y, and z and returns the
// maximum cosine distance.
//
// cos distance = 1.0 - a.dot(b) / (norm(a) * norm(b))
double ComputeCosineDistance(const Eigen::MatrixXd& mat1,
                             const Eigen::MatrixXd& mat2) {
  Eigen::Vector3d cos_distance;
  for (int i =0; i < 3; i++) {
    cos_distance(i) = 1.0 - std::abs(mat1.row(i).dot(mat2.row(i)));
  }
  return cos_distance.maxCoeff();
}

// Find the maximal rigid component containing fixed_node. This is done by
// examining which nodes are parallel when removing node fixed_node from the
// null space. The nodes are only parallel if they are part of the maximal rigid
// component with fixed_node.
void FindMaximalParallelRigidComponent(const Eigen::MatrixXd& null_space,
                                       const int fixed_node,
                                       std::unordered_set<int>* largest_cc) {
  static const double kMaxCosDistance = 1e-5;
  static const double kMaxNorm = 1e-10;

  const int num_nodes = null_space.rows() / 3;

  largest_cc->insert(fixed_node);

  const Eigen::MatrixXd fixed_null_space_component =
    null_space.block(3 * fixed_node, 0, 3, null_space.cols());

  // Remove the fixed node from the rest of the null space.
  Eigen::MatrixXd modified_null_space =
      null_space - fixed_null_space_component.replicate(num_nodes, 1);

  // Normalize all rows to be unit-norm. If the rows have a very small norm then
  // they are parallel to the fixed_node. and should be set to zero.
  const Eigen::VectorXd norms = modified_null_space.rowwise().norm();

  modified_null_space.rowwise().normalize();

  // Find the pairs to match. Add all indices that are close to 0-vectors, as
  // they are clearly part of the rigid component.
  std::vector<int> indices_to_match;
  for (int i = 0; i < num_nodes; i++) {
    // Skip this index if it is fixed.
    if (i == fixed_node) {
      continue;
    }

    // Skip this index if it is nearly a 0-vector because this means it is
    // clearly part of the rigid component.
    if (norms(3 * i) < kMaxNorm &&
        norms(3 * i + 1) < kMaxNorm &&
        norms(3 * i + 2) < kMaxNorm) {
      largest_cc->insert(i);
      continue;
    }

    indices_to_match.emplace_back(i);
  }

  // Each node has three dimensions (x, y, z). We only compare parallel-ness
  // between similar dimensions. If all x, y, z dimensions are parallel then
  // the two nodes will be parallel.
  for (int i = 0; i < indices_to_match.size(); i++) {
    // Test all other nodes (that have not been tested) to determine if they are
    // parallel to this node.
    const Eigen::MatrixXd& block1 = modified_null_space.block(
        3 * indices_to_match[i], 0, 3, null_space.cols());

    for (int j = i + 1; j < indices_to_match.size(); j++) {
      const Eigen::MatrixXd& block2 = modified_null_space.block(
          3 * indices_to_match[j], 0, 3, null_space.cols());

      const double cos_distance = ComputeCosineDistance(block1, block2);
      if (cos_distance < kMaxCosDistance) {
        largest_cc->insert(indices_to_match[i]);
        largest_cc->insert(indices_to_match[j]);
      }
    }
  }
}

}  // namespace

void ExtractMaximallyParallelRigidSubgraph(
    const std::unordered_map<ViewId, Eigen::Vector3d>& orientations,
    ViewGraph* view_graph) {
  // Create a mapping of indexes to ViewIds for our linear system.
  std::unordered_map<ViewId, int> view_ids_to_index;
  view_ids_to_index.reserve(orientations.size());
  for (const auto& orientation : orientations) {
    if (!view_graph->HasView(orientation.first)) {
      continue;
    }
    const int current_index = view_ids_to_index.size();
    InsertIfNotPresent(&view_ids_to_index, orientation.first, current_index);
  }

  // Form the global angle measurements matrix from:
  //    t_{i,j} x (c_j - c_i) = 0.
  Eigen::MatrixXd angle_measurements(3 * view_graph->NumEdges(),
                                     3 * orientations.size());
  FormAngleMeasurementMatrix(orientations,
                             *view_graph,
                             view_ids_to_index,
                             &angle_measurements);

  // Extract the null space of the angle measurements matrix.
  Eigen::FullPivLU<Eigen::MatrixXd> lu(angle_measurements.transpose() *
                                       angle_measurements);
  const Eigen::MatrixXd null_space = lu.kernel();

  // For each node in the graph (i.e. each camera), set the null space component
  // to be zero such that the camera position would be fixed at the origin. If
  // two nodes i and j are in the same rigid component, then their null spaces
  // will be parallel because the camera positions may only change by a
  // scale. We find all components that are parallel to find the rigid
  // components. The largest of such component is the maximally parallel rigid
  // component of the graph.
  std::unordered_set<int> maximal_rigid_component;
  for (int i = 0; i < orientations.size(); i++) {
    std::unordered_set<int> temp_cc;
    FindMaximalParallelRigidComponent(null_space, i, &temp_cc);
    if (temp_cc.size() > maximal_rigid_component.size()) {
      std::swap(temp_cc, maximal_rigid_component);
    }
  }

  // Only keep the nodes in the largest maximally parallel rigid component.
  for (const auto& orientation : orientations) {
    const int index = FindOrDie(view_ids_to_index, orientation.first);
    // If the view is not in the maximal rigid component then remove it from the
    // view graph.
    if (!ContainsKey(maximal_rigid_component, index) &&
        view_graph->HasView(orientation.first)) {
      CHECK(view_graph->RemoveView(orientation.first))
          << "Could not remove view id " << orientation.first
          << " from the view graph because it does not exist.";
    }
  }
}

}  // namespace theia
