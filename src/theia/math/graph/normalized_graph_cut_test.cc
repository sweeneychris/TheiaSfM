// Copyright (C) 2015 The Regents of the University of California (Regents).
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

#include <unordered_map>
#include <unordered_set>
#include <utility>

#include "gtest/gtest.h"
#include "theia/math/graph/normalized_graph_cut.h"
#include "theia/util/hash.h"

namespace theia {

// Create a simple scenario with 4 nodes that form a rectangle:
//    0 ------------------------- 3
//    |                           |
//    1 ------------------------- 4
//    |                           |
//    2 ------------------------- 5
// This should be very simple to partition.
TEST(NormalizedGraphCut, SimpleGraph) {
//  typedef std::pair<int, int> IntPair;
//  std::unordered_map<std::pair<int, int>, double> edge_weights;
//  edge_weights.emplace(IntPair(0, 1), 1000);
//  edge_weights.emplace(IntPair(1, 2), 1000);
//  edge_weights.emplace(IntPair(0, 3), 1000);
//  edge_weights.emplace(IntPair(3, 4), 1000);
//  edge_weights.emplace(IntPair(4, 5), 1000);
//  edge_weights.emplace(IntPair(3, 5), 1000);
//  edge_weights.emplace(IntPair(0, 3), 1);
//  edge_weights.emplace(IntPair(1, 4), 1);
//  edge_weights.emplace(IntPair(2, 5), 1);
//
//  NormalizedGraphCut<int>::Options options;
//  NormalizedGraphCut<int> ncut(options);
//  std::unordered_set<int> subgraph1, subgraph2;
//  EXPECT_TRUE(ncut.ComputeCut(edge_weights, &subgraph1, &subgraph2, NULL));
//
//  // Make sure that the subgraphs are split properly.
//  EXPECT_EQ(subgraph1.size(), 3);
//  EXPECT_EQ(subgraph2.size(), 3);
//
//  const int node_0_subgraph = ContainsKey(subgraph1, 0) ? 1 : 2;
//  const int node_1_subgraph = ContainsKey(subgraph1, 1) ? 1 : 2;
//  const int node_2_subgraph = ContainsKey(subgraph1, 2) ? 1 : 2;
//  const int node_3_subgraph = ContainsKey(subgraph1, 3) ? 1 : 2;
//  const int node_4_subgraph = ContainsKey(subgraph1, 3) ? 1 : 2;
//  const int node_5_subgraph = ContainsKey(subgraph1, 3) ? 1 : 2;
//
//  EXPECT_EQ(node_0_subgraph, node_1_subgraph);
//  EXPECT_EQ(node_1_subgraph, node_2_subgraph);
//  EXPECT_EQ(node_3_subgraph, node_4_subgraph);
//  EXPECT_EQ(node_4_subgraph, node_5_subgraph);
//  EXPECT_NE(node_0_subgraph, node_3_subgraph);
}

}  // namespace theia
