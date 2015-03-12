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

#include "theia/sfm/filter_view_pairs_from_orientation.h"

#include <ceres/rotation.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <glog/logging.h>
#include <unordered_map>
#include <unordered_set>

#include "theia/math/util.h"
#include "theia/util/hash.h"
#include "theia/util/map_util.h"
#include "theia/sfm/twoview_info.h"
#include "theia/sfm/types.h"
#include "theia/sfm/view_graph/view_graph.h"

namespace theia {

namespace {

bool AngularDifferenceIsAcceptable(
    const Eigen::Vector3d& orientation1,
    const Eigen::Vector3d& orientation2,
    const Eigen::Vector3d& relative_orientation,
    const double max_relative_rotation_difference_degrees) {
  Eigen::Matrix3d rotation_matrix1, rotation_matrix2, relative_rotation_matrix;
  ceres::AngleAxisToRotationMatrix(
      orientation1.data(),
      ceres::ColumnMajorAdapter3x3(rotation_matrix1.data()));
  ceres::AngleAxisToRotationMatrix(
      orientation2.data(),
      ceres::ColumnMajorAdapter3x3(rotation_matrix2.data()));
  ceres::AngleAxisToRotationMatrix(
      relative_orientation.data(),
      ceres::ColumnMajorAdapter3x3(relative_rotation_matrix.data()));

  // Compose the relative rotation from the two known orientation.
  const Eigen::Matrix3d composed_relative_rotation_matrix =
      rotation_matrix2 * rotation_matrix1.transpose();

  // Compute angular distance between the relative rotations.
  const double rotation_angular_difference_degrees =
      RadToDeg(Eigen::AngleAxisd(relative_rotation_matrix.transpose() *
                                 composed_relative_rotation_matrix).angle());
  return rotation_angular_difference_degrees <=
         max_relative_rotation_difference_degrees;
}

}  // namespace

void FilterViewPairsFromOrientation(
    const std::unordered_map<ViewId, Eigen::Vector3d>& orientations,
    const double max_relative_rotation_difference_degrees,
    ViewGraph* view_graph) {
  CHECK_NOTNULL(view_graph);
  CHECK_GE(max_relative_rotation_difference_degrees, 0.0);

  std::unordered_set<ViewIdPair> view_pairs_to_remove;
  const auto& view_pairs = view_graph->GetAllEdges();
  for (const auto& view_pair : view_pairs) {
    const Eigen::Vector3d* orientation1 =
        FindOrNull(orientations, view_pair.first.first);
    const Eigen::Vector3d* orientation2 =
        FindOrNull(orientations, view_pair.first.second);

    // If the view pair contains a view that does not have an orientation then
    // remove it.
    if (orientation1 == nullptr || orientation2 == nullptr) {
      LOG(WARNING)
          << "View pair (" << view_pair.first.first << ", "
          << view_pair.first.second
          << ") contains a view that does not exist! Removing the view pair.";
      view_pairs_to_remove.insert(view_pair.first);
      continue;
    }

    // Remove the view pair if the relative rotation estimate is not within the
    // tolerance.
    if (!AngularDifferenceIsAcceptable(
            *orientation1,
            *orientation2,
            view_pair.second.rotation_2,
            max_relative_rotation_difference_degrees)) {
      view_pairs_to_remove.insert(view_pair.first);
    }
  }

  // Remove all the "bad" relative poses.
  for (const ViewIdPair view_id_pair : view_pairs_to_remove) {
    view_graph->RemoveEdge(view_id_pair.first, view_id_pair.second);
  }
  VLOG(1) << "Removed " << view_pairs_to_remove.size()
          << " view pairs by rotation filtering.";
}

}  // namespace theia
