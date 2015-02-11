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

#include "theia/sfm/view_graph/triplet_extractor.h"

#include <glog/logging.h>
#include <algorithm>
#include <unordered_map>
#include <memory>
#include <utility>
#include <vector>

#include "theia/math/graph/connected_components.h"
#include "theia/sfm/types.h"
#include "theia/sfm/view_triplet.h"

namespace theia {

// Extracts all triplets from the view pairs. Triplets are grouped by
// connectivity, and vector represents a connected triplet graph.
bool TripletExtractor::ExtractTripletsFromViewPairs(
    const std::unordered_map<ViewIdPair, TwoViewInfo>& view_pairs,
    std::vector<std::vector<ViewTriplet> >* connected_triplets) {
  triplets_.clear();
  triplet_edges_.clear();
  triplets_.reserve(view_pairs.size());
  triplet_edges_.reserve(view_pairs.size());

  // Find all the triplets.
  FindTripletsInViewPairs(view_pairs);

  // Add the triplets to the connected component.
  ConnectedComponents<TripletId> triplets_cc_extractor;
  AddTripletsToConnectedComponents(&triplets_cc_extractor);

  // Extract the connected components.
  std::unordered_map<TripletId, std::unordered_set<TripletId> >
      triplets_connected_components;
  triplets_cc_extractor.Extract(&triplets_connected_components);

  // Move the connected triplets to the output.
  connected_triplets->reserve(triplets_connected_components.size());
  for (const auto& connected_component : triplets_connected_components) {
    VLOG(2) << "Extracted a connected triplet graph of containing "
            << connected_component.second.size() << " triplet(s)";

    std::vector<ViewTriplet> triplets;
    triplets.reserve(connected_component.second.size());
    for (const auto& triplet_id : connected_component.second) {
      triplets.emplace_back(triplets_[triplet_id]);
    }
    connected_triplets->emplace_back(triplets);
  }
  return true;
}

// Finds all triplets in the view pairs. We find loops by examining a sorted
// list of the view ids of the edges in the view graph.
void TripletExtractor::FindTripletsInViewPairs(
    const std::unordered_map<ViewIdPair, TwoViewInfo>& view_pairs) {
  // Get a sorted list of view pair ids.
  std::vector<ViewIdPair> view_pair_ids;
  view_pair_ids.reserve(view_pairs.size());
  for (const auto& view_pair : view_pairs) {
    view_pair_ids.emplace_back(view_pair.first);
  }
  std::sort(view_pair_ids.begin(), view_pair_ids.end());

  // Iterate the list to find triplets. For a given vertex, we choose all
  // combinations of edges that contain that vertex and check if the connecting
  // vertices contain an edge between them. For example, given the sorted view
  // pair ids:
  //
  //  {0, 1}
  //  {0, 3}
  //  {0, 4}
  //  {1, 2}
  //  {1, 3}
  //   ...
  //
  // we would take a pair of edges containing vertex 0, such as {0, 1} and {0,
  // 3} then check to see if the edge {1, 3} exists. If done on a sorted list,
  // this should be very fast.

  for (int i = 0; i < view_pair_ids.size(); i++) {
    for (int j = i + 1;
         j < view_pair_ids.size() &&
             view_pair_ids[i].first == view_pair_ids[j].first;
         j++) {
      // A triplet only exists if a connecting edge exists.
      const ViewIdPair edge_1_3_id(view_pair_ids[i].second,
                                   view_pair_ids[j].second);
      if (ContainsKey(view_pairs, edge_1_3_id)) {
        const ViewId triplet_view_ids[3] = {view_pair_ids[i].first,
                                            view_pair_ids[i].second,
                                            view_pair_ids[j].second};
        TwoViewInfo triplet_edges[3];
        triplet_edges[0] = FindOrDieNoPrint(view_pairs, view_pair_ids[i]);
        triplet_edges[1] = FindOrDieNoPrint(view_pairs, view_pair_ids[j]);
        triplet_edges[2] = FindOrDieNoPrint(view_pairs, edge_1_3_id);
        AddViewTriplet(triplet_view_ids, triplet_edges);
      }
    }
  }
}

// Creates and adds a view triplet to the triplets_ container.
void TripletExtractor::AddViewTriplet(const ViewId view_ids[3],
                                      const TwoViewInfo edges[3]) {
  // Add triplet to the triplet container.
  ViewTriplet triplet;
  triplet.view_ids[0] = view_ids[0];
  triplet.view_ids[1] = view_ids[1];
  triplet.view_ids[2] = view_ids[2];
  triplet.info_one_two = edges[0];
  triplet.info_one_three = edges[1];
  triplet.info_two_three = edges[2];
  triplets_.emplace_back(triplet);

  const TripletId triplet_id = triplets_.size() - 1;

  // Add each edge of the triplet to the lookup map.
  triplet_edges_[ViewIdPair(view_ids[0], view_ids[1])].insert(triplet_id);
  triplet_edges_[ViewIdPair(view_ids[0], view_ids[2])].insert(triplet_id);
  triplet_edges_[ViewIdPair(view_ids[1], view_ids[2])].insert(triplet_id);
}

// Each view triplet contains 3 view pairs, so we use the lookup map to add
// all triplets that share one of the view pairs as neighbors in the connected
// component analysis.
void TripletExtractor::AddTripletsToConnectedComponents(
    ConnectedComponents<TripletId>* triplet_cc) {
  for (const auto& triplet_edge : triplet_edges_) {
    // Add all possible triplets that share this view pair.
    for (auto iterator1 = triplet_edge.second.begin();
         iterator1 != triplet_edge.second.end();
         ++iterator1) {
      // We add one edge with the triplet connected to itself so that singleton
      // triplets are still considered.
      for (auto iterator2 = iterator1; iterator2 != triplet_edge.second.end();
           ++iterator2) {
        triplet_cc->AddEdge(*iterator1, *iterator2);
      }
    }
  }
}

}  // namespace theia
