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

#ifndef THEIA_SFM_VIEW_GRAPH_TRIPLET_EXTRACTOR_H_
#define THEIA_SFM_VIEW_GRAPH_TRIPLET_EXTRACTOR_H_

#include <stdint.h>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "theia/util/util.h"
#include "theia/util/hash.h"
#include "theia/sfm/types.h"
#include "theia/sfm/view_triplet.h"

namespace theia {
class TwoViewInfo;
template<typename T> class ConnectedComponents;

// Extract all loops of size 3 (i.e., triplets) in a set of view pairs. Triplets
// are then gathered into connected components where two triplets are connected
// if the share an edge in the view pairs. NOTE: This means that a single
// connected view graph may results in multiple "connected" triplet graphs.
class TripletExtractor {
 public:
  typedef uint32_t TripletId;

  TripletExtractor() {}

  // Extracts all triplets from the view pairs (which should be edges in a view
  // graph). Triplets are grouped by connectivity, and vector represents a
  // connected triplet graph.
  bool ExtractTripletsFromViewPairs(
      const std::unordered_map<ViewIdPair, TwoViewInfo>& view_pairs,
      std::vector<std::vector<ViewTriplet> >* connected_triplets);

 private:
  // Finds all triplets in the view pairs.
  void FindTripletsInViewPairs(
      const std::unordered_map<ViewIdPair, TwoViewInfo>& view_graph);

  // Creates and adds a view triplet to the triplets_ container.
  void AddViewTriplet(const ViewId view_ids[3], const TwoViewInfo edges[3]);

  // Each view triplet contains 3 view pairs, so we use the lookup map to add
  // all triplets that share one of the view pairs as neighbors in the connected
  // component analysis.
  void AddTripletsToConnectedComponents(
      ConnectedComponents<TripletId>* triplet_cc);

  // Container for all triplets found in the view pairs.
  std::vector<ViewTriplet> triplets_;

  // We keep track of the triplets that each edge in the view pairs participate
  // in for computing the connected components.
  std::unordered_map<ViewIdPair, std::unordered_set<TripletId> > triplet_edges_;

  DISALLOW_COPY_AND_ASSIGN(TripletExtractor);
};

}  // namespace theia

#endif  // THEIA_SFM_VIEW_GRAPH_TRIPLET_EXTRACTOR_H_
