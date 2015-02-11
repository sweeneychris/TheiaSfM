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

#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "gtest/gtest.h"
#include "theia/math/graph/connected_components.h"

namespace theia {

// Fully connected graph.
TEST(ConnectedComponents, FullyConnectedGraph) {
  ConnectedComponents<int> connected_components;

  for (int i = 0; i < 9; i++) {
    connected_components.AddEdge(i, i + 1);
  }

  std::unordered_map<int, std::unordered_set<int> > disjoint_sets;
  connected_components.Extract(&disjoint_sets);
  EXPECT_EQ(disjoint_sets.size(), 1);
  EXPECT_EQ(disjoint_sets.begin()->second.size(), 10);
}

// Fully disconnected graph.
TEST(ConnectedComponents, FullyDisconnectedGraph) {
  ConnectedComponents<int> connected_components;

  for (int i = 0; i < 9; i++) {
    connected_components.AddEdge(i, i + 100);
  }

  std::unordered_map<int, std::unordered_set<int> > disjoint_sets;
  connected_components.Extract(&disjoint_sets);
  EXPECT_EQ(disjoint_sets.size(), 9);
  for (const auto& cc : disjoint_sets) {
    EXPECT_EQ(cc.second.size(), 2);
  }
}

// Partially connected graph.
TEST(ConnectedComponents, PartiallyConnectedGraph) {
  ConnectedComponents<int> connected_components;

  for (int i = 0; i < 9; i++) {
    connected_components.AddEdge(i, i + 5);
  }

  std::unordered_map<int, std::unordered_set<int> > disjoint_sets;
  connected_components.Extract(&disjoint_sets);
  EXPECT_EQ(disjoint_sets.size(), 5);
  for (const auto& cc : disjoint_sets) {
    EXPECT_GE(cc.second.size(), 2);
  }
}

// Graph with size limitation.
TEST(ConnectedComponents, LimitComponentSize) {
  ConnectedComponents<int> connected_components(2);

  for (int i = 0; i < 9; i++) {
    connected_components.AddEdge(i, i + 1);
  }

  std::unordered_map<int, std::unordered_set<int> > disjoint_sets;
  connected_components.Extract(&disjoint_sets);
  EXPECT_EQ(disjoint_sets.size(), 5);
  for (const auto& cc : disjoint_sets) {
    EXPECT_LE(cc.second.size(), 2);
  }
}

// Random shuffle in order does not matter.
TEST(ConnectedComponents, RandomOrder) {
  std::vector<std::pair<int, int> > pairs_to_add;
  for (int i = 0; i < 10; i++) {
    pairs_to_add.emplace_back(i, i + 1);
  }

  for (int i = 0; i < 25; i++) {
    ConnectedComponents<int> connected_components;
    std::random_shuffle(pairs_to_add.begin(), pairs_to_add.end());

    for (const auto& pair_to_add : pairs_to_add) {
      connected_components.AddEdge(pair_to_add.first, pair_to_add.second);
    }

    std::unordered_map<int, std::unordered_set<int> > disjoint_sets;
    connected_components.Extract(&disjoint_sets);
    EXPECT_EQ(disjoint_sets.size(), 1);
  }
}


}  // namespace theia
