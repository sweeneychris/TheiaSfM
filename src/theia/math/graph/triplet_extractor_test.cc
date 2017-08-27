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
#include <tuple>
#include <utility>
#include <vector>

#include "gtest/gtest.h"

#include "theia/math/graph/triplet_extractor.h"

namespace theia {
typedef std::pair<int, int> IntPair;
typedef std::tuple<int, int, int> IntTriplet;

TEST(ViewTriplet, NoTriplets) {
  std::unordered_set<IntPair> edges;
  edges.emplace(0, 1);
  edges.emplace(1, 2);
  edges.emplace(2, 3);

  TripletExtractor<int> triplet_extractor;
  std::vector<std::vector<IntTriplet> > triplets;
  triplet_extractor.ExtractTriplets(edges, &triplets);
  EXPECT_EQ(triplets.size(), 0);
}

TEST(ViewTriplet, TwoTripletsOneSet) {
  std::unordered_set<IntPair> edges;
  edges.emplace(0, 1);
  edges.emplace(1, 2);
  edges.emplace(2, 3);
  edges.emplace(0, 3);
  edges.emplace(1, 3);
  // Add an edge to the view graph that will not be part of any triplets.
  edges.emplace(3, 4);

  TripletExtractor<int> triplet_extractor;
  std::vector<std::vector<IntTriplet> > triplets;
  triplet_extractor.ExtractTriplets(edges, &triplets);

  EXPECT_EQ(triplets.size(), 1);
  EXPECT_EQ(triplets.at(0).size(), 2);
}

TEST(ViewTriplet, DisconnectedSets) {
  std::unordered_set<IntPair> edges;
  edges.emplace(0, 1);
  edges.emplace(1, 2);
  edges.emplace(0, 2);
  edges.emplace(0, 3);
  edges.emplace(0, 4);
  edges.emplace(3, 4);

  TripletExtractor<int> triplet_extractor;
  std::vector<std::vector<IntTriplet> > triplets;
  triplet_extractor.ExtractTriplets(edges, &triplets);

  EXPECT_EQ(triplets.size(), 2);
  EXPECT_EQ(triplets.at(0).size(), 1);
  EXPECT_EQ(triplets.at(1).size(), 1);
}

}  // namespace theia
