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

#include "theia/sfm/view_graph/orientations_from_maximum_spanning_tree.h"

#include <Eigen/Core>
#include <unordered_map>
#include <unordered_set>

#include "theia/math/graph/minimum_spanning_tree.h"
#include "theia/sfm/global_pose_estimation/linear_rotation_estimator.h"
#include "theia/sfm/twoview_info.h"
#include "theia/sfm/types.h"
#include "theia/sfm/view_graph/view_graph.h"
#include "theia/util/map_util.h"

namespace theia {

bool OrientationsFromMaximumSpanningTree(
    const ViewGraph& view_graph,
    std::unordered_map<ViewId, Eigen::Vector3d>* orientations) {
  CHECK_NOTNULL(orientations);

  // Compute maximum spanning tree.
  const auto& all_edges = view_graph.GetAllEdges();
  MinimumSpanningTree<ViewId, int> mst_extractor;
  for (const auto& edge : all_edges) {
    // Since we want the *maximum* spanning tree, we negate all of the edge
    // weights in the *minimum* spanning tree extractor.
    mst_extractor.AddEdge(edge.first.first,
                          edge.first.second,
                          -edge.second.num_verified_matches);
  }

  std::unordered_set<ViewIdPair> mst;
  if (!mst_extractor.Extract(&mst)) {
    VLOG(2)
        << "Could not extract the maximum spanning tree from the view graph";
    return false;
  }

  // Collect the MST edges for the rotation computation.
  std::unordered_map<ViewIdPair, TwoViewInfo> mst_view_pairs;
  mst_view_pairs.reserve(mst.size());
  for (const auto& edge : mst) {
    mst_view_pairs[edge] = FindOrDieNoPrint(all_edges, edge);
  }

  // Solve for rotations using the linear solver. Weigh the edges by the inlier
  // count to give more weight to more stable edges.
  LinearRotationEstimator rotation_estimator(true);
  return rotation_estimator.EstimateRotations(mst_view_pairs, orientations);
}

}  // namespace theia
