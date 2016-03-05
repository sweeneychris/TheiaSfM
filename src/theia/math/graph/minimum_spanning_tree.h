// Copyright (C) 2016 The Regents of the University of California (Regents).
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

#ifndef THEIA_MATH_GRAPH_MINIMUM_SPANNING_TREE_H_
#define THEIA_MATH_GRAPH_MINIMUM_SPANNING_TREE_H_

#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "theia/math/graph/connected_components.h"
#include "theia/util/hash.h"
#include "theia/util/util.h"

namespace theia {

// A class for extracting the minimum spanning tree of a graph using Kruskal's
// greedy algorithm. The minimum spanning tree is a subgraph that contains all
// nodes in the graph and only the edges that connect these nodes with a minimum
// edge weight summation. The algorithm runs in O(E * log (V) ) where E is the
// number of edges and V is the number of nodes in the graph. For more details
// on the algorithm please see:
//   https://en.wikipedia.org/wiki/Kruskal%27s_algorithm
template <typename T, typename V>
class MinimumSpanningTree {
 public:
  MinimumSpanningTree() {}

  // Add an edge in the graph.
  void AddEdge(const T& node1, const T& node2, const V& weight) {
    edges_.emplace_back(weight, std::pair<T, T>(node1, node2));
  }

  // Extracts the minimum spanning tree. Returns true on success and false upon
  // failure. If true is returned, the output variable contains the edge list of
  // the minimum spanning tree.
  bool Extract(std::unordered_set<std::pair<T, T> >* minimum_spanning_tree) {
    if (edges_.size() == 0) {
      VLOG(2) << "No edges were passed to the minimum spanning tree extractor!";
      return false;
    }

    // Determine the number of nodes in the graph.
    const int num_nodes = CountNodesInGraph();

    // Reserve space in the MST since we know it will have exactly N - 1 edges.
    minimum_spanning_tree->reserve(num_nodes - 1);

    // Order all edges by their weights.
    std::sort(edges_.begin(), edges_.end());

    // For each edge in the graph, add it to the minimum spanning tree if it
    // does not create a cycle.
    ConnectedComponents<T> cc;
    for (int i = 0;
         i < edges_.size() && minimum_spanning_tree->size() < num_nodes - 1;
         i++) {
      const auto& edge = edges_[i];
      if (!cc.NodesInSameConnectedComponent(edge.second.first,
                                            edge.second.second)) {
        cc.AddEdge(edge.second.first, edge.second.second);
        minimum_spanning_tree->emplace(edge.second.first, edge.second.second);
      }
    }

    return minimum_spanning_tree->size() == num_nodes - 1;
  }

 private:
  // Counts the number of nodes in the graph by counting the number of unique
  // node values we have received from AddEdge.
  int CountNodesInGraph() {
    std::vector<T> nodes;
    nodes.reserve(edges_.size() * 2);
    for (const auto& edge : edges_) {
      nodes.emplace_back(edge.second.first);
      nodes.emplace_back(edge.second.second);
    }
    std::sort(nodes.begin(), nodes.end());
    auto unique_end = std::unique(nodes.begin(), nodes.end());
    return std::distance(nodes.begin(), unique_end);
  }

  std::vector<std::pair<V, std::pair<T, T> > > edges_;

  // Each node is mapped to a Root node. If the node is equal to the root id
  // then the node is a root and the size of the root is the size of the
  // connected component.
  std::unordered_map<T, T> disjoint_set_;

  DISALLOW_COPY_AND_ASSIGN(MinimumSpanningTree);
};

}  // namespace theia

#endif  // THEIA_MATH_GRAPH_MINIMUM_SPANNING_TREE_H_
