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

#include <memory>
#include <vector>
#include "gtest/gtest.h"

#include "theia/sfm/view_graph/triplet_extractor.h"
#include "theia/sfm/view_graph/view_graph.h"
#include "theia/sfm/view_triplet.h"

namespace theia {
namespace {

// Add a triplet (composed of 3 view pairs) to the view graph. The two view info
// is not used so we set it to the default.
void AddTripletToViewGraph(const ViewId view_id1,
                           const ViewId view_id2,
                           const ViewId view_id3,
                           ViewGraph* view_graph) {
  TwoViewInfo info;
  view_graph->AddEdge(view_id1, view_id2, info);
  view_graph->AddEdge(view_id1, view_id3, info);
  view_graph->AddEdge(view_id2, view_id3, info);
}

}  // namespace

TEST(ViewTriplet, NoTriplets) {
  ViewGraph view_graph;
  TwoViewInfo info;
  view_graph.AddEdge(0, 1, info);
  view_graph.AddEdge(1, 2, info);
  view_graph.AddEdge(2, 3, info);

  const auto& view_pairs = view_graph.GetAllEdges();
  TripletExtractor triplet_extractor;
  std::vector<std::vector<ViewTriplet> > triplets;
  triplet_extractor.ExtractTripletsFromViewPairs(view_pairs, &triplets);
  EXPECT_EQ(triplets.size(), 0);
}

TEST(ViewTriplet, TwoTripletsOneSet) {
  ViewGraph view_graph;
  AddTripletToViewGraph(0, 1, 2, &view_graph);
  AddTripletToViewGraph(0, 1, 3, &view_graph);

  // Add an edge to the view graph that will not be part of any triplets.
  TwoViewInfo info;
  view_graph.AddEdge(3, 4, info);

  const auto& view_pairs = view_graph.GetAllEdges();
  TripletExtractor triplet_extractor;
  std::vector<std::vector<ViewTriplet> > triplets;
  triplet_extractor.ExtractTripletsFromViewPairs(view_pairs, &triplets);

  EXPECT_EQ(triplets.size(), 1);
  EXPECT_EQ(triplets.at(0).size(), 2);
}

TEST(ViewTriplet, DisconnectedSets) {
  ViewGraph view_graph;
  AddTripletToViewGraph(0, 1, 2, &view_graph);
  AddTripletToViewGraph(0, 3, 4, &view_graph);

  const auto& view_pairs = view_graph.GetAllEdges();
  TripletExtractor triplet_extractor;
  std::vector<std::vector<ViewTriplet> > triplets;
  triplet_extractor.ExtractTripletsFromViewPairs(view_pairs, &triplets);

  EXPECT_EQ(triplets.size(), 2);
  EXPECT_EQ(triplets.at(0).size(), 1);
  EXPECT_EQ(triplets.at(1).size(), 1);
}

}  // namespace theia
