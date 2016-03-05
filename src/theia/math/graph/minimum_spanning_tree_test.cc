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

#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "gtest/gtest.h"
#include "theia/math/graph/connected_components.h"
#include "theia/math/graph/minimum_spanning_tree.h"
#include "theia/util/map_util.h"
#include "theia/util/random.h"

namespace theia {

typedef std::pair<int, int> IntPair;

void VerifyMST(const std::unordered_map<IntPair, double>& edges,
               const std::unordered_set<IntPair>& mst,
               const double mst_cost) {
  // Check that all nodes were added using connected component.
  std::unordered_set<int> nodes;
  for (const auto& edge : edges) {
    nodes.emplace(edge.first.first);
    nodes.emplace(edge.first.second);
  }

  ConnectedComponents<int> cc_mst;
  for (const auto& edge : mst) {
    cc_mst.AddEdge(edge.first, edge.second);
  }
  std::unordered_map<int, std::unordered_set<int> > components_mst;
  cc_mst.Extract(&components_mst);

  // Verify that there is only 1 connected component in the MST.
  EXPECT_EQ(components_mst.size(), 1);

  // Verify that the MST covers all nodes.
  const auto& cc_iterator = components_mst.begin();
  EXPECT_EQ(nodes.size(), cc_iterator->second.size());
  for (const auto& node : cc_iterator->second) {
    EXPECT_TRUE(ContainsKey(nodes, node));
  }

  // Verify it is the MST.
  double cost = 0;
  for (const auto& edge : mst) {
    cost += FindOrDieNoPrint(edges, edge);
  }

  EXPECT_EQ(cost, mst_cost);
}

TEST(MinimumSpanningTree, SimpleSpanningGraph) {
  MinimumSpanningTree<int, double> mst;
  std::unordered_map<IntPair, double> edges;
  for (int i = 0; i < 9; i++) {
    mst.AddEdge(i, i + 1, 1.0);
    edges.emplace(IntPair(i, i + 1), 1.0);
  }
  const double mst_cost = edges.size() * 1.0;

  std::unordered_set<IntPair> spanning_tree;
  EXPECT_TRUE(mst.Extract(&spanning_tree));
  VerifyMST(edges, spanning_tree, mst_cost);
}

TEST(MinimumSpanningTree, RandomGraph) {
  const int num_vertices = 10;
  const int num_edges = 30;

  MinimumSpanningTree<int, double> mst;
  std::unordered_map<IntPair, double> edges;

  // Create the MST first.
  for (int i = 0; i < num_vertices - 1; i++) {
    mst.AddEdge(i, i + 1, 1.0);
    edges.emplace(IntPair(i, i + 1), 1.0);
  }
  const double mst_cost = edges.size() * 1.0;

  // Add extra edges with a higher edge weight.
  while (edges.size() != num_edges) {
    const int random_node1 = RandInt(0, num_vertices - 1);
    const int random_node2 = RandInt(0, num_vertices - 1);
    if (random_node1 == random_node2 ||
        ContainsKey(edges, IntPair(random_node1, random_node2)) ||
        ContainsKey(edges, IntPair(random_node2, random_node1))) {
      continue;
    }

    const double edge_weight = RandDouble(2.0, 10.0);
    mst.AddEdge(random_node1, random_node2, edge_weight);
    edges.emplace(IntPair(random_node1, random_node2), edge_weight);
  }

  std::unordered_set<IntPair> spanning_tree;
  EXPECT_TRUE(mst.Extract(&spanning_tree));
  VerifyMST(edges, spanning_tree, mst_cost);
}

TEST(MinimumSpanningTree, NoEdges) {
  MinimumSpanningTree<int, double> mst;
  std::unordered_set<IntPair> spanning_tree;
  EXPECT_FALSE(mst.Extract(&spanning_tree));
}

TEST(MinimumSpanningTree, DisconnectedGraph) {
  // Create a graph: 0 -> 1 -> 2  and 3 -> 4 -> 5
  MinimumSpanningTree<int, double> mst;
  mst.AddEdge(0, 1, 1.0);
  mst.AddEdge(1, 2, 1.0);
  mst.AddEdge(3, 4, 1.0);
  mst.AddEdge(4, 5, 1.0);
  std::unordered_set<IntPair> spanning_tree;
  EXPECT_FALSE(mst.Extract(&spanning_tree));
}

}  // namespace theia
