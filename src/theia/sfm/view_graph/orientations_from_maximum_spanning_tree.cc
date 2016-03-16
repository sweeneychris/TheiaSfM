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

#include <ceres/rotation.h>
#include <Eigen/Core>

#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "theia/math/graph/minimum_spanning_tree.h"
#include "theia/sfm/twoview_info.h"
#include "theia/sfm/types.h"
#include "theia/sfm/view_graph/view_graph.h"
#include "theia/util/map_util.h"

namespace theia {
namespace {
typedef std::pair<TwoViewInfo, ViewIdPair> HeapElement;

bool SortHeapElement(const HeapElement& h1, const HeapElement& h2) {
  return h1.first.num_verified_matches > h2.first.num_verified_matches;
}

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

// Adds all the edges of view_id to the heap. Only edges that do not already
// have an orientation estimation are added.
void AddEdgesToHeap(
    const ViewGraph& view_graph,
    const std::unordered_map<ViewId, Eigen::Vector3d>& orientations,
    const ViewId view_id,
    std::vector<HeapElement >* heap) {
  const std::unordered_set<ViewId>* edge_ids =
      view_graph.GetNeighborIdsForView(view_id);
  for (const ViewId edge_id : *edge_ids) {
    // Only add edges to the heap that contain a vertex that has not been seen.
    if (ContainsKey(orientations, edge_id)) {
      continue;
    }

    heap->emplace_back(*view_graph.GetEdge(view_id, edge_id),
                       ViewIdPair(view_id, edge_id));
    std::push_heap(heap->begin(), heap->end(), SortHeapElement);
  }
}

}  // namespace

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

  // Create an MST view graph.
  ViewGraph mst_view_graph;
  for (const ViewIdPair& edge : mst) {
    mst_view_graph.AddEdge(edge.first, edge.second,
                           *view_graph.GetEdge(edge.first, edge.second));
  }

  // Chain the relative rotations together to compute orientations.  We use a
  // heap to determine the next edges to add to the minimum spanning tree.
  std::vector<HeapElement> heap;

  // Set the root value.
  const ViewId root_view_id = mst.begin()->first;
  (*orientations)[root_view_id] = Eigen::Vector3d::Zero();
  AddEdgesToHeap(mst_view_graph, *orientations, root_view_id, &heap);

  while (!heap.empty()) {
    const HeapElement next_edge = heap.front();
    // Remove the best edge.
    std::pop_heap(heap.begin(), heap.end(), SortHeapElement);
    heap.pop_back();

    // If the edge contains two vertices that have already been added then do
    // nothing.
    if (ContainsKey(*orientations, next_edge.second.second)) {
      continue;
    }

    // Compute the orientation for the vertex.
    (*orientations)[next_edge.second.second] = ComputeOrientation(
        FindOrDie(*orientations, next_edge.second.first),
        next_edge.first,
        next_edge.second.first,
        next_edge.second.second);

    // Add all edges to the heap.
    AddEdgesToHeap(mst_view_graph,
                   *orientations,
                   next_edge.second.second,
                   &heap);
  }
  return true;
}

}  // namespace theia
