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

#include "theia/sfm/view_graph/orientations_from_view_graph.h"

#include <ceres/rotation.h>
#include <Eigen/Core>
#include <queue>
#include <unordered_map>

#include "theia/util/map_util.h"
#include "theia/sfm/types.h"
#include "theia/sfm/view_graph/view_graph.h"

namespace theia {

namespace {

// Computes the orientation of the neighbor camera based on the orientation of
// the source camera and the relative rotation between the cameras.
Eigen::Vector3d ComputeOrientation(const Eigen::Vector3d& source_orientation,
                                   const TwoViewInfo& two_view_info,
                                   const ViewId source_view_id,
                                   const ViewId neighbor_view_id) {
  Eigen::Matrix3d source_rotation_mat, relative_rotation;
  ceres::AngleAxisToRotationMatrix(
      source_orientation.data(),
      ceres::ColumnMajorAdapter3x3(source_rotation_mat.data()));
  ceres::AngleAxisToRotationMatrix(
      two_view_info.rotation_2.data(),
      ceres::ColumnMajorAdapter3x3(relative_rotation.data()));

  const Eigen::Matrix3d neighbor_orientation =
      (source_view_id < neighbor_view_id)
          ? (relative_rotation * source_rotation_mat).eval()
          : (relative_rotation.transpose() * source_rotation_mat).eval();

  Eigen::Vector3d orientation;
  ceres::RotationMatrixToAngleAxis(
      ceres::ColumnMajorAdapter3x3(neighbor_orientation.data()),
      orientation.data());
  return orientation;
}

}  // namespace

void OrientationsFromViewGraph(
    const ViewGraph& view_graph,
    const ViewId root_view_id,
    std::unordered_map<ViewId, Eigen::Vector3d>* orientations) {
  CHECK(view_graph.HasView(root_view_id))
      << "The root node does not exist in the view graph.";

  std::queue<ViewId> views_to_visit;

  // Set the root value and push it to the queue.
  (*orientations)[root_view_id] = Eigen::Vector3d::Zero();
  views_to_visit.push(root_view_id);
  while (!views_to_visit.empty()) {
    const ViewId view_id = views_to_visit.front();
    views_to_visit.pop();
    const Eigen::Vector3d& source_orientation =
        FindOrDie(*orientations, view_id);

    // Add all the neigbors that have not been visited.
    const auto* neighbors = view_graph.GetEdgesForView(view_id);
    for (const auto& neighbor : *neighbors) {
      if (ContainsKey(*orientations, neighbor.first)) {
        continue;
      }

      // Compute the orientation and add it to the result.
      (*orientations)[neighbor.first] = ComputeOrientation(source_orientation,
                                                           neighbor.second,
                                                           view_id,
                                                           neighbor.first);
      // Add the neighbor to the queue.
      views_to_visit.push(neighbor.first);
    }
  }
}

}  // namespace theia
