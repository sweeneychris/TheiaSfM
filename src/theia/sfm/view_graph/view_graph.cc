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

#include <cereal/archives/portable_binary.hpp>
#include <cstdio>
#include <cstdlib>
#include <fstream>   // NOLINT
#include <iostream>  // NOLINT
#include <unordered_map>
#include <unordered_set>

#include "theia/sfm/twoview_info.h"
#include "theia/sfm/types.h"
#include "theia/util/hash.h"
#include "theia/util/map_util.h"

namespace theia {

// Number of views in the graph.
int ViewGraph::NumViews() const { return vertices_.size(); }

int ViewGraph::NumEdges() const { return edges_.size(); }

// Utilities to read and write a view graph to/from disk.
bool ViewGraph::ReadFromDisk(const std::string& input_file) {
  std::ifstream input_reader(input_file, std::ios::in | std::ios::binary);
  if (!input_reader.is_open()) {
    LOG(ERROR) << "Could not open the file: " << input_file << " for reading.";
    return false;
  }

  cereal::PortableBinaryInputArchive input_archive(input_reader);
  input_archive(*this);

  return true;
}

bool ViewGraph::WriteToDisk(const std::string& output_file) {
  std::ofstream output_writer(output_file, std::ios::out | std::ios::binary);
  if (!output_writer.is_open()) {
    LOG(ERROR) << "Could not open the file: " << output_file << " for writing.";
    return false;
  }

  cereal::PortableBinaryOutputArchive output_archive(output_writer);
  output_archive(*this);

  return true;
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

bool ViewGraph::HasEdge(const ViewId view_id_1, const ViewId view_id_2) const {
  const ViewIdPair view_id_pair = (view_id_1 < view_id_2)
                                      ? ViewIdPair(view_id_1, view_id_2)
                                      : ViewIdPair(view_id_2, view_id_1);
  return ContainsKey(edges_, view_id_pair);
}

// Removes the view from the view graph and removes all edges connected to the
// view. Returns true on success and false if the view did not exist in the
// view graph.
bool ViewGraph::RemoveView(const ViewId view_id) {
  const auto* neighbor_ids = FindOrNull(vertices_, view_id);
  if (neighbor_ids == nullptr) {
    return false;
  }

  // Remove the edges to the view from adjacent vertices.
  for (const ViewId neighbor_id : *neighbor_ids) {
    vertices_[neighbor_id].erase(view_id);
    const ViewIdPair view_id_pair = (view_id < neighbor_id)
                                        ? ViewIdPair(view_id, neighbor_id)
                                        : ViewIdPair(neighbor_id, view_id);
    edges_.erase(view_id_pair);
  }

  // Remove the view as a vertex.
  vertices_.erase(view_id);
  return true;
}

// Adds an edge between the two views with the edge value of
// two_view_info. New vertices are added to the graph if they did not already
// exist. If an edge already existed between the two views then the edge value
// is updated.
void ViewGraph::AddEdge(const ViewId view_id_1,
                        const ViewId view_id_2,
                        const TwoViewInfo& two_view_info) {
  if (view_id_1 == view_id_2) {
    DLOG(WARNING) << "Cannot add an edge from view id " << view_id_1
                  << " to itself!";
    return;
  }

  const ViewIdPair view_id_pair = (view_id_1 < view_id_2)
                                      ? ViewIdPair(view_id_1, view_id_2)
                                      : ViewIdPair(view_id_2, view_id_1);

  DLOG_IF(WARNING, ContainsKey(edges_, view_id_pair))
      << "An edge already exists between view " << view_id_1 << " and view "
      << view_id_2;

  vertices_[view_id_1].insert(view_id_2);
  vertices_[view_id_2].insert(view_id_1);
  edges_[view_id_pair] = two_view_info;
}

// Removes the edge from the view graph. Returns true if the edge is removed
// and false if the edge did not exist.
bool ViewGraph::RemoveEdge(const ViewId view_id_1, const ViewId view_id_2) {
  const ViewIdPair view_id_pair = (view_id_1 < view_id_2)
                                      ? ViewIdPair(view_id_1, view_id_2)
                                      : ViewIdPair(view_id_2, view_id_1);
  if (!ContainsKey(edges_, view_id_pair)) {
    return false;
  }

  // Erase the edge from each vertex.
  if (vertices_[view_id_1].erase(view_id_2) == 0 ||
      vertices_[view_id_2].erase(view_id_1) == 0 ||
      edges_.erase(view_id_pair) == 0) {
    return false;
  }

  return true;
}

// Returns all the edges for a given
const std::unordered_set<ViewId>* ViewGraph::GetNeighborIdsForView(
    const ViewId view_id) const {
  return FindOrNull(vertices_, view_id);
}

// Returns the edge value or NULL if it does not exist.
const TwoViewInfo* ViewGraph::GetEdge(const ViewId view_id_1,
                                      const ViewId view_id_2) const {
  const ViewIdPair view_id_pair = (view_id_1 < view_id_2)
                                      ? ViewIdPair(view_id_1, view_id_2)
                                      : ViewIdPair(view_id_2, view_id_1);
  return FindOrNull(edges_, view_id_pair);
}

TwoViewInfo* ViewGraph::GetMutableEdge(const ViewId view_id_1,
                                       const ViewId view_id_2) {
  const ViewIdPair view_id_pair = (view_id_1 < view_id_2)
                                      ? ViewIdPair(view_id_1, view_id_2)
                                      : ViewIdPair(view_id_2, view_id_1);
  return FindOrNull(edges_, view_id_pair);
}

// Returns a map of all edges. Each edge is found exactly once in the map and
// is indexed by the ViewIdPair (view id 1, view id 2) such that view id 1 <
// view id 2.
const std::unordered_map<ViewIdPair, TwoViewInfo>& ViewGraph::GetAllEdges()
    const {
  return edges_;
}

// Extract a subgraph containing only the specified views.
void ViewGraph::ExtractSubgraph(
    const std::unordered_set<ViewId>& views_in_subgraph,
    ViewGraph* subgraph) const {
  CHECK_NOTNULL(subgraph);

  // Iterate over each vertex and add all edges to that vertex that are part of
  // the subgraph.
  for (const auto& vertex : vertices_) {
    // If the vertex is not contained in the subgraph, skip it.
    if (!ContainsKey(views_in_subgraph, vertex.first)) {
      continue;
    }

    // Iterate over all edges to the current vertex and find edges to add to the
    // subgraph. An edge will be added to the subgraph only if both vertices are
    // in the subgraph.
    for (const ViewId& second_vertex : vertex.second) {
      // Skip this vertex (and thus, edge) if it is not in the subgraph. Also,
      // only consider edges where vertex1 < vertex2 so as not to add redundant
      // to the subgraph.
      if (!ContainsKey(views_in_subgraph, second_vertex) ||
          second_vertex < vertex.first) {
        continue;
      }

      // Add the edge to the subgraph.
      const TwoViewInfo& edge =
          FindOrDieNoPrint(edges_, ViewIdPair(vertex.first, second_vertex));
      subgraph->AddEdge(vertex.first, second_vertex, edge);
    }
  }
}

}  // namespace theia
