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
#include <algorithm>
#include <unordered_map>

#include "theia/util/map_util.h"
#include "theia/sfm/twoview_info.h"
#include "theia/sfm/types.h"
#include "theia/sfm/view_graph/view_graph.h"

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
  const auto* edges = view_graph.GetEdgesForView(view_id);
  for (const auto& edge : *edges) {
    // Only add edges to the heap that contain a vertex that has not been seen.
    if (ContainsKey(orientations, edge.first)) {
      continue;
    }

    heap->emplace_back(edge.second, ViewIdPair(view_id, edge.first));
    std::push_heap(heap->begin(), heap->end(), SortHeapElement);
  }
}

}  // namespace

void OrientationsFromViewGraph(
    const ViewGraph& view_graph,
    const ViewId root_view_id,
    std::unordered_map<ViewId, Eigen::Vector3d>* orientations) {
  CHECK(view_graph.HasView(root_view_id))
      << "The root node does not exist in the view graph.";

  // We use a heap to determine the next edges to add to the minimum spanning
  // tree.
  std::vector<HeapElement> heap;

  // Set the root value.
  (*orientations)[root_view_id] = Eigen::Vector3d::Zero();
  AddEdgesToHeap(view_graph, *orientations, root_view_id, &heap);

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
    AddEdgesToHeap(view_graph, *orientations, next_edge.second.second, &heap);
  }
}

}  // namespace theia
