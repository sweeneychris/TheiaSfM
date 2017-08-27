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

#ifndef THEIA_MATH_GRAPH_TRIPLET_EXTRACTOR_H_
#define THEIA_MATH_GRAPH_TRIPLET_EXTRACTOR_H_

#include <stdint.h>

#include <algorithm>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "theia/math/graph/connected_components.h"
#include "theia/sfm/types.h"
#include "theia/util/hash.h"
#include "theia/util/util.h"

namespace theia {
// Extract all loops of size 3 (i.e., triplets) in a set of view pairs. Triplets
// are then gathered into connected components where two triplets are connected
// if the share an edge in the view pairs. NOTE: This means that a single
// connected view graph may results in multiple "connected" triplet graphs.
template <typename T>
class TripletExtractor {
 public:
  typedef uint32_t TripletId;
  typedef std::pair<T, T> TypePair;
  typedef std::tuple<T, T, T> TypeTriplet;

  TripletExtractor() {}

  // Extracts all triplets from the view pairs (which should be edges in a view
  // graph). Triplets are grouped by connectivity, and vector represents a
  // connected triplet graph.
  bool ExtractTriplets(
      const std::unordered_set<TypePair>& edge_graph,
      std::vector<std::vector<TypeTriplet>>* connected_triplets);

 private:
  // Finds all triplets in the view pairs.
  void FindTriplets(const std::unordered_set<TypePair>& edge_graph);

  // Each view triplet contains 3 view pairs, so we use the lookup map to add
  // all triplets that share one of the view pairs as neighbors in the connected
  // component analysis.
  void GetConnectedTripletGraphs(
      std::unordered_map<TripletId, std::unordered_set<TripletId>>*
          connected_triplet_graphs);

  // Store the triplet in the internal container and add the entries to the edge
  // list appropriatesly.
  void StoreTriplet(const T& a, const T& b, const T& c);

  // Container for all triplets found in the view pairs.
  std::vector<TypeTriplet> triplets_;

  // We keep track of the triplets that each edge in the view pairs participate
  // in for computing the connected components.
  std::unordered_map<TypePair, std::vector<TripletId>> triplet_edges_;

  DISALLOW_COPY_AND_ASSIGN(TripletExtractor);
};

namespace internal {
// Compare the two containers by size such that the larger container moves to
// the front of a sorted vector.
template <typename T>
bool CompareBySize(const std::vector<std::tuple<T, T, T>>& v1,
                   const std::vector<std::tuple<T, T, T>>& v2) {
  return v1.size() > v2.size();
}

// The following logic will sort the triplet in increasing order.
template <typename T>
void SortTriplet(std::tuple<T, T, T>* tuple) {
  T& a = std::get<0>(*tuple);
  T& b = std::get<1>(*tuple);
  T& c = std::get<2>(*tuple);
  // We know that a < b (by construction) so we only need to correctly place c.
  if (c < b) {
    std::swap(b, c);
  }
  if (b < a) {
    std::swap(a, b);
  }
}
}  // namespace internal

template <typename T>
bool TripletExtractor<T>::ExtractTriplets(
    const std::unordered_set<TypePair>& edge_graph,
    std::vector<std::vector<TypeTriplet>>* connected_triplets) {
  triplets_.clear();
  triplet_edges_.clear();
  triplets_.reserve(edge_graph.size());
  triplet_edges_.reserve(edge_graph.size());

  // Find all the triplets.
  FindTriplets(edge_graph);

  // Split the triplets into connected triplet graphs.
  std::unordered_map<TripletId, std::unordered_set<TripletId>>
      connected_triplet_graphs;
  GetConnectedTripletGraphs(&connected_triplet_graphs);

  // Move the connected triplets to the output.
  connected_triplets->reserve(connected_triplet_graphs.size());
  for (const auto& connected_component : connected_triplet_graphs) {
    VLOG(2) << "Extracted a connected triplet graph of containing "
            << connected_component.second.size() << " triplet(s)";

    std::vector<TypeTriplet> triplets;
    triplets.reserve(connected_component.second.size());
    for (const auto& triplet_id : connected_component.second) {
      triplets.emplace_back(triplets_[triplet_id]);
    }
    connected_triplets->emplace_back(triplets);
  }

  // Sort the triplet connected components such that the largest is at the front
  // of the output.
  std::sort(connected_triplets->begin(),
            connected_triplets->end(),
            internal::CompareBySize<T>);

  return true;
}

// Finds all triplets in the view pairs. We find loops by examining a sorted
// list of the view ids of the edges in the view graph.
template <typename T>
void TripletExtractor<T>::FindTriplets(
    const std::unordered_set<TypePair>& edge_graph) {
  // Create an adjacency graph.
  std::unordered_map<T, std::set<T>> adjacency_graph;
  adjacency_graph.reserve(edge_graph.size());
  for (const TypePair& edge : edge_graph) {
    adjacency_graph[edge.first].emplace(edge.second);
    adjacency_graph[edge.second].emplace(edge.first);
  }

  // Iterate over the adjacency graph and find all triplets.
  for (const auto& edge : edge_graph) {
    // Find any nodes that contain edges to both of the nodes for the current
    // edge. Any such nodes form a triplet in the graph. We do this efficiently
    // using a set intersection on the list of adjacent nodes to each of the
    // nodes in the current edges.
    const std::set<T>& edges_to_node1 =
        theia::FindOrDie(adjacency_graph, edge.first);
    const std::set<T>& edges_to_node2 =
        theia::FindOrDie(adjacency_graph, edge.second);

    // Compute the intersection between the two adjacency lists to find
    // triplets.
    std::vector<T> node_intersection;
    node_intersection.reserve(
        std::min(edges_to_node1.size(), edges_to_node2.size()));
    std::set_intersection(edges_to_node1.begin(),
                          edges_to_node1.end(),
                          edges_to_node2.begin(),
                          edges_to_node2.end(),
                          std::back_inserter(node_intersection));
    // Keep track of each triplet.
    for (int i = 0; i < node_intersection.size(); i++) {
      StoreTriplet(edge.first, edge.second, node_intersection[i]);
    }
    // Remove the edge from the adjacency graph. Since this step found all
    // triplets that contain the edge, there are no possible triplets remaining
    // that include the edge so it may be safely removed to improve efficiency.
    adjacency_graph[edge.first].erase(edge.second);
    adjacency_graph[edge.second].erase(edge.first);
  }
}

template <typename T>
void TripletExtractor<T>::StoreTriplet(const T& a, const T& b, const T& c) {
  TypeTriplet triplet(a, b, c);
  internal::SortTriplet(&triplet);
  triplets_.emplace_back(triplet);

  // Add each edge of the triplet to the edge lookup map. This will be used
  // to determine connected components.
  const TripletId triplet_id = triplets_.size() - 1;

  // Add each edge of the triplet to the lookup map.
  const TypePair pair_01(std::get<0>(triplet), std::get<1>(triplet));
  const TypePair pair_02(std::get<0>(triplet), std::get<2>(triplet));
  const TypePair pair_12(std::get<1>(triplet), std::get<2>(triplet));
  triplet_edges_[pair_01].emplace_back(triplet_id);
  triplet_edges_[pair_02].emplace_back(triplet_id);
  triplet_edges_[pair_12].emplace_back(triplet_id);
}

// Each view triplet contains 3 view pairs, so we use the lookup map to add
// all triplets that share one of the view pairs as neighbors in the connected
// component analysis.
template <typename T>
void TripletExtractor<T>::GetConnectedTripletGraphs(
    std::unordered_map<TripletId, std::unordered_set<TripletId>>*
        connected_triplet_graphs) {
  ConnectedComponents<TripletId> triplets_cc_extractor;

  // Iterate through all edges and connect the triplets that contain the edge.
  for (const auto& triplet_edge : triplet_edges_) {
    const std::vector<TripletId>& triplets_containing_edge =
        triplet_edge.second;
    // We add a self-contained edge to a triplet if it is in isolation. This
    // ensures that the triplet is still returned in the CC analysis even though
    // it has no neighboring triplets.
    if (triplets_containing_edge.size() == 1) {
      triplets_cc_extractor.AddEdge(triplets_containing_edge[0],
                                    triplets_containing_edge[0]);
      continue;
    }
    // For connected component analysis, we simply need to link all triplets
    // containing this edge together so they will be in the same conneted
    // component. We can do this by simply chaining them together.
    for (int i = 0; i < triplets_containing_edge.size() - 1; i++) {
      triplets_cc_extractor.AddEdge(triplets_containing_edge[i],
                                    triplets_containing_edge[i + 1]);
    }
  }

  triplets_cc_extractor.Extract(connected_triplet_graphs);
}

}  // namespace theia

#endif  // THEIA_MATH_GRAPH_TRIPLET_EXTRACTOR_H_
