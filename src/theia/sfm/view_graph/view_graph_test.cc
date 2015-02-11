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

#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "gtest/gtest.h"

#include "theia/util/hash.h"
#include "theia/util/map_util.h"
#include "theia/sfm/twoview_info.h"
#include "theia/sfm/types.h"
#include "theia/sfm/view_graph/view_graph.h"

namespace theia {

// This is needed for EXPECT_EQ(TwoViewInfo, TwoViewInfo);
bool operator==(const TwoViewInfo& lhs, const TwoViewInfo& rhs) {
  return lhs.position_2 == rhs.position_2 &&
      lhs.rotation_2 == rhs.rotation_2 &&
      lhs.num_verified_matches == rhs.num_verified_matches;
}

TEST(ViewGraph, Constructor) {
  ViewGraph view_graph;
  EXPECT_EQ(view_graph.NumViews(), 0);
  EXPECT_EQ(view_graph.NumEdges(), 0);
}


TEST(ViewGraph, AddView) {
  ViewGraph graph;
  const ViewId view_one = 0;
  const ViewId view_two = 1;

  graph.AddView(view_one);
  EXPECT_TRUE(graph.HasView(view_one));
  EXPECT_EQ(graph.NumViews(), 1);
  EXPECT_EQ(graph.NumEdges(), 0);

  graph.AddView(view_two);
  EXPECT_TRUE(graph.HasView(view_two));
  EXPECT_EQ(graph.NumViews(), 2);
  EXPECT_EQ(graph.NumEdges(), 0);

  graph.AddView(view_one);
  EXPECT_EQ(graph.NumViews(), 2);
  EXPECT_EQ(graph.NumEdges(), 0);
}

TEST(ViewGraph, RemoveView) {
  ViewGraph graph;
  const ViewId view_one = 0;
  const ViewId view_two = 1;
  const TwoViewInfo edge;

  graph.AddView(view_one);
  graph.AddView(view_two);
  graph.AddEdge(view_one, view_two, edge);
  EXPECT_TRUE(graph.HasView(view_one));
  EXPECT_TRUE(graph.HasView(view_two));
  EXPECT_EQ(graph.NumEdges(), 1);

  graph.RemoveView(view_one);
  EXPECT_TRUE(!graph.HasView(view_one));
  EXPECT_EQ(graph.NumViews(), 1);
  EXPECT_EQ(graph.NumEdges(), 0);
}

TEST(ViewGraph, AddEdge) {
  TwoViewInfo info1;
  info1.num_verified_matches = 1;
  TwoViewInfo info2;
  info2.num_verified_matches = 2;

  ViewGraph graph;
  graph.AddEdge(0, 1, info1);
  graph.AddEdge(0, 2, info2);

  EXPECT_EQ(graph.NumViews(), 3);
  EXPECT_EQ(graph.NumEdges(), 2);

  const std::vector<int> expected_ids = {0, 1, 2};
  std::unordered_set<ViewId> view_ids = graph.ViewIds();

  const auto* edge_0_1 = graph.GetEdge(0, 1);
  const auto* edge_0_2 = graph.GetEdge(0, 2);
  EXPECT_TRUE(edge_0_1 != nullptr);
  EXPECT_EQ(*edge_0_1, info1);
  EXPECT_TRUE(edge_0_2 != nullptr);
  EXPECT_EQ(*edge_0_2, info2);
  EXPECT_TRUE(graph.GetEdge(1, 2) == nullptr);

  const auto* edge_1_0 = graph.GetEdge(1, 0);
  const auto* edge_2_0 = graph.GetEdge(2, 0);
  EXPECT_TRUE(edge_1_0 != nullptr);
  EXPECT_EQ(*edge_1_0, info1);
  EXPECT_EQ(*edge_2_0, info2);
  EXPECT_TRUE(graph.GetEdge(2, 1) == nullptr);

  const std::unordered_map<ViewId, TwoViewInfo> edges =
      *graph.GetEdgesForView(0);
  EXPECT_TRUE(ContainsKey(edges, 1));
  EXPECT_EQ(FindOrDie(edges, 1), info1);
  EXPECT_TRUE(ContainsKey(edges, 2));
  EXPECT_EQ(FindOrDie(edges, 2), info2);

  TwoViewInfo info3;
  info3.num_verified_matches = 3;
  graph.AddEdge(1, 2, info3);
  EXPECT_EQ(graph.NumViews(), 3);
  EXPECT_EQ(graph.NumEdges(), 3);
  EXPECT_TRUE(graph.GetEdge(1, 2) != nullptr);
  EXPECT_EQ(*graph.GetEdge(1, 2), info3);
}

TEST(ViewGraph, RemoveEdge) {
  TwoViewInfo info1;
  info1.num_verified_matches = 1;
  TwoViewInfo info2;
  info2.num_verified_matches = 2;

  ViewGraph graph;
  graph.AddEdge(0, 1, info1);
  graph.AddEdge(0, 2, info2);

  EXPECT_TRUE(graph.GetEdge(0, 1) != nullptr);
  EXPECT_TRUE(graph.GetEdge(0, 2) != nullptr);

  EXPECT_EQ(graph.NumViews(), 3);
  EXPECT_EQ(graph.NumEdges(), 2);

  EXPECT_TRUE(graph.RemoveEdge(0, 2));
  EXPECT_EQ(graph.NumEdges(), 1);
  EXPECT_TRUE(graph.GetEdge(0, 2) == nullptr);
}

}  // namespace theia
