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

#include "theia/sfm/view_graph/view_graph.h"

#include <unordered_map>
#include <unordered_set>

#include "theia/util/hash.h"
#include "theia/util/map_util.h"
#include "theia/sfm/twoview_info.h"
#include "theia/sfm/types.h"

namespace theia {

// Number of views in the graph.
int ViewGraph::NumViews() const { return vertices_.size(); }

int ViewGraph::NumEdges() const {
  int num_edges = 0;
  for (const auto& vertex : vertices_) {
    num_edges += vertex.second.size();
  }
  // Divide by 2 because each undirected edge will be counted twice.
  return num_edges / 2;
}

// Returns a set of the ViewIds contained in the view graph.
std::unordered_set<ViewId> ViewGraph::ViewIds() const {
  std::unordered_set<ViewId> view_ids;
  view_ids.reserve(vertices_.size());
  for (const auto& vertex : vertices_) {
    view_ids.insert(vertex.first);
  }
  return view_ids;
}

bool ViewGraph::HasView(const ViewId view_id) const {
  return ContainsKey(vertices_, view_id);
}

// Adds the view to the view graph if it is not already present.
void ViewGraph::AddView(const ViewId view_id) {
  if (ContainsKey(vertices_, view_id)) {
    DLOG(INFO) << "View " << view_id << " already exists in view graph!";
    return;
  }
  vertices_[view_id] = std::unordered_map<ViewId, TwoViewInfo>();
}

// Removes the view from the view graph and removes all edges connected to the
// view. Returns true on success and false if the view did not exist in the
// view graph.
bool ViewGraph::RemoveView(const ViewId view_id) {
  const auto* edges = FindOrNull(vertices_, view_id);
  if (edges == nullptr) {
    return false;
  }

  // Remove the edges to the view from adjacent vertices.
  for (const auto& adjacent_view : *edges) {
    vertices_[adjacent_view.first].erase(view_id);
  }

  // Remove the view as a vertex.
  vertices_.erase(view_id);
  return true;
}

// Adds an edge between the two views with the edge value of
// two_view_info. New vertices are added to the graph if they did not already
// exist. If an edge already existed between the two views then the edge value
// is updated.
void ViewGraph::AddEdge(const ViewId view_id_1, const ViewId view_id_2,
                        const TwoViewInfo& two_view_info) {
  if (view_id_1 == view_id_2) {
    DLOG(WARNING) << "Cannot add an edge from view id " << view_id_1
                  << " to itself!";
    return;
  }

  DLOG_IF(WARNING, ContainsKey(vertices_[view_id_1], view_id_2))
      << "An edge already exists between view " << view_id_1 << " and view "
      << view_id_2;

  vertices_[view_id_1][view_id_2] = two_view_info;
  vertices_[view_id_2][view_id_1] = two_view_info;
}

// Removes the edge from the view graph. Returns true if the edge is removed
// and false if the edge did not exist.
bool ViewGraph::RemoveEdge(const ViewId view_id_1, const ViewId view_id_2) {
  if (!ContainsKey(vertices_, view_id_1) ||
      !ContainsKey(vertices_, view_id_2)) {
    return false;
  }

  // Erase the edge from each vertex.
  if (vertices_[view_id_1].erase(view_id_2) == 0 ||
      vertices_[view_id_2].erase(view_id_1) == 0) {
    return false;
  }

  return true;
}

// Returns all the edges for a given
const std::unordered_map<ViewId, TwoViewInfo>* ViewGraph::GetEdgesForView(
    const ViewId view_id) const {
  return FindOrNull(vertices_, view_id);
}

// Returns the edge value or NULL if it does not exist.
const TwoViewInfo* ViewGraph::GetEdge(const ViewId view_id_1,
                                      const ViewId view_id_2) const {
  auto* edges = FindOrNull(vertices_, view_id_1);
  return (edges != nullptr) ? FindOrNull(*edges, view_id_2) : nullptr;
}

// Returns a map of all edges. Each edge is found exactly once in the map and
// is indexed by the ViewIdPair (view id 1, view id 2) such that view id 1 <
// view id 2.
std::unordered_map<ViewIdPair, TwoViewInfo> ViewGraph::GetAllEdges() const {
  std::unordered_map<ViewIdPair, TwoViewInfo> two_view_infos;
  two_view_infos.reserve(NumEdges());
  for (const auto& vertex : vertices_) {
    for (const auto& adjacent_vertex : vertex.second) {
      // Skip the edge if view 1 > view 2. This prevents us from visiting the
      // same edge twice.
      if (vertex.first > adjacent_vertex.first) {
        continue;
      }

      two_view_infos[ViewIdPair(vertex.first, adjacent_vertex.first)] =
          adjacent_vertex.second;
    }
  }
  return two_view_infos;
}

}  // namespace theia
