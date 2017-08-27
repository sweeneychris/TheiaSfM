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

#include "theia/sfm/filter_view_graph_cycles_by_rotation.h"

#include <Eigen/Core>
#include <ceres/rotation.h>
#include <glog/logging.h>
#include <unordered_map>
#include <vector>

#include "theia/math/graph/triplet_extractor.h"
#include "theia/math/util.h"
#include "theia/sfm/twoview_info.h"
#include "theia/sfm/types.h"
#include "theia/sfm/view_graph/view_graph.h"
#include "theia/sfm/view_triplet.h"
#include "theia/util/hash.h"

namespace theia {
namespace {

double ComputeLoopRotationError(const TwoViewInfo& info_one_two,
                                const TwoViewInfo& info_one_three,
                                const TwoViewInfo& info_two_three) {
  // Get relative rotation matrices.
  Eigen::Matrix3d rotation1_2, rotation1_3, rotation2_3;
  ceres::AngleAxisToRotationMatrix(
      info_one_two.rotation_2.data(),
      ceres::ColumnMajorAdapter3x3(rotation1_2.data()));
  ceres::AngleAxisToRotationMatrix(
      info_one_three.rotation_2.data(),
      ceres::ColumnMajorAdapter3x3(rotation1_3.data()));
  ceres::AngleAxisToRotationMatrix(
      info_two_three.rotation_2.data(),
      ceres::ColumnMajorAdapter3x3(rotation2_3.data()));

  // Compute loop rotation.
  const Eigen::Matrix3d loop_rotation =
      rotation2_3 * rotation1_2 * rotation1_3.transpose();
  Eigen::Vector3d loop_rotation_angle_axis;
  ceres::RotationMatrixToAngleAxis(
      ceres::ColumnMajorAdapter3x3(loop_rotation.data()),
      loop_rotation_angle_axis.data());

  // Return the angle of the loop rotation which is the error of the
  // concatenated triplet rotation.
  return RadToDeg(loop_rotation_angle_axis.norm());
}

void RemoveTripletEdgesFromInvalidViewPairs(
    const ViewIdTriplet& triplet,
    std::unordered_set<ViewIdPair>* invalid_view_pairs) {
  const ViewIdPair pair_01(std::get<0>(triplet), std::get<1>(triplet));
  const ViewIdPair pair_02(std::get<0>(triplet), std::get<2>(triplet));
  const ViewIdPair pair_12(std::get<1>(triplet), std::get<2>(triplet));
  invalid_view_pairs->erase(pair_01);
  invalid_view_pairs->erase(pair_02);
  invalid_view_pairs->erase(pair_12);
}

}  // namespace

void FilterViewGraphCyclesByRotation(const double max_loop_error_degrees,
                                     ViewGraph* view_graph) {
  const std::unordered_map<ViewIdPair, TwoViewInfo>& view_pairs =
      view_graph->GetAllEdges();

  // Initialize a list of invalid view pairs to all view pairs. View pairs
  // deemed valid will be removed from this list.
  std::unordered_set<ViewIdPair> invalid_view_pairs;
  invalid_view_pairs.reserve(view_pairs.size());
  for (const auto& view_pair : view_pairs) {
    invalid_view_pairs.insert(view_pair.first);
  }

  // Find all triplets.
  TripletExtractor<ViewId> triplet_extractor;
  std::vector<std::vector<ViewIdTriplet> > connected_triplets;
  CHECK(triplet_extractor.ExtractTriplets(invalid_view_pairs,
                                          &connected_triplets))
      << "Could not extract triplets from view pairs.";

  // Examine the cycles of size 3 to determine invalid view pairs from the
  // rotations.
  for (const auto& triplets : connected_triplets) {
    for (const ViewIdTriplet& triplet : triplets) {
      const TwoViewInfo& info_one_two =
          *view_graph->GetEdge(std::get<0>(triplet), std::get<1>(triplet));
      const TwoViewInfo& info_one_three =
          *view_graph->GetEdge(std::get<0>(triplet), std::get<2>(triplet));
      const TwoViewInfo& info_two_three =
          *view_graph->GetEdge(std::get<1>(triplet), std::get<2>(triplet));

      // Compute loop rotation error.
      const double loop_rotation_error_degrees = ComputeLoopRotationError(
          info_one_two, info_one_three, info_two_three);

      // Add the view pairs to the list of valid view pairs if the loop error is
      // within the designated tolerance.
      if (loop_rotation_error_degrees < max_loop_error_degrees) {
        RemoveTripletEdgesFromInvalidViewPairs(triplet, &invalid_view_pairs);
      }
    }
  }

  VLOG(1) << "Removing " << invalid_view_pairs.size() << " of "
          << view_pairs.size() << " view pairs from loop rotation filtering.";

  // Remove any view pairs not in the list of valid edges.
  for (const ViewIdPair& view_pair_id : invalid_view_pairs) {
    view_graph->RemoveEdge(view_pair_id.first, view_pair_id.second);
  }
}

}  // namespace theia
