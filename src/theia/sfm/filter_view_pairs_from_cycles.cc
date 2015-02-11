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

#include "theia/sfm/filter_view_pairs_from_cycles.h"

#include <ceres/rotation.h>
#include <Eigen/Core>
#include <glog/logging.h>
#include <unordered_map>
#include <vector>

#include "theia/math/util.h"
#include "theia/util/hash.h"
#include "theia/sfm/twoview_info.h"
#include "theia/sfm/types.h"
#include "theia/sfm/view_graph/triplet_extractor.h"
#include "theia/sfm/view_triplet.h"

namespace theia {

namespace {

double ComputeLoopRotationError(const ViewTriplet& triplet) {
  // Get relative rotation matrices.
  Eigen::Matrix3d rotation1_2, rotation1_3, rotation2_3;
  ceres::AngleAxisToRotationMatrix(
      triplet.info_one_two.rotation_2.data(),
      ceres::ColumnMajorAdapter3x3(rotation1_2.data()));
  ceres::AngleAxisToRotationMatrix(
      triplet.info_one_three.rotation_2.data(),
      ceres::ColumnMajorAdapter3x3(rotation1_3.data()));
  ceres::AngleAxisToRotationMatrix(
      triplet.info_two_three.rotation_2.data(),
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
    const ViewTriplet& triplet,
    std::unordered_set<ViewIdPair>* invalid_view_pairs) {
  invalid_view_pairs->erase(
      ViewIdPair(triplet.view_ids[0], triplet.view_ids[1]));
  invalid_view_pairs->erase(
      ViewIdPair(triplet.view_ids[0], triplet.view_ids[2]));
  invalid_view_pairs->erase(
      ViewIdPair(triplet.view_ids[1], triplet.view_ids[2]));
}

}  // namespace

void FilterViewPairsFromCycles(
    const double max_loop_error_degrees,
    std::unordered_map<ViewIdPair, TwoViewInfo>* view_pairs) {
  // Find all triplets.
  TripletExtractor triplet_extractor;
  std::vector<std::vector<ViewTriplet> > connected_triplets;
  CHECK(triplet_extractor.ExtractTripletsFromViewPairs(*view_pairs,
                                                       &connected_triplets))
      << "Could not extract triplets from view pairs.";

  // Create a list of invalid view pairs. View pairs deemed valid will be
  // removed from this list.
  std::unordered_set<ViewIdPair> invalid_view_pairs;
  invalid_view_pairs.reserve(view_pairs->size());
  for (const auto& view_pair : *view_pairs) {
    invalid_view_pairs.insert(view_pair.first);
  }

  // Examine the cycles of size 3 to determine invalid view pairs from the
  // rotations.
  for (const auto& triplets : connected_triplets) {
    for (const ViewTriplet& triplet : triplets) {
      // Compute loop rotation error.
      const double loop_rotation_error_degrees =
          ComputeLoopRotationError(triplet);

      // Add the view pairs to the list of valid view pairs if the loop error is
      // within the designated tolerance.
      if (loop_rotation_error_degrees < max_loop_error_degrees) {
        RemoveTripletEdgesFromInvalidViewPairs(triplet, &invalid_view_pairs);
      }
    }
  }

  VLOG(1) << "Removing " << invalid_view_pairs.size() << " of "
          << view_pairs->size() << " view pairs from loop rotation filtering.";

  // Remove any view pairs not in the list of valid edges.
  for (const ViewIdPair& view_pair_id : invalid_view_pairs) {
    view_pairs->erase(view_pair_id);
  }
}

}  // namespace theia
