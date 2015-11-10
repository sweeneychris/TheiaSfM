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

#include "theia/sfm/view_graph/remove_disconnected_view_pairs.h"

#include <glog/logging.h>
#include <unordered_map>
#include <unordered_set>

#include "theia/math/graph/connected_components.h"
#include "theia/sfm/types.h"
#include "theia/sfm/view_graph/view_graph.h"

namespace theia {

// Removes all view pairs that are not part of the largest connected component.
std::unordered_set<ViewId> RemoveDisconnectedViewPairs(ViewGraph* view_graph) {
  CHECK_NOTNULL(view_graph);
  std::unordered_set<ViewId> removed_views;

  // Extract all connected components.
  ConnectedComponents<ViewId> cc_extractor;
  const auto& view_pairs = view_graph->GetAllEdges();
  for (const auto& view_pair : view_pairs) {
    cc_extractor.AddEdge(view_pair.first.first, view_pair.first.second);
  }
  std::unordered_map<ViewId, std::unordered_set<ViewId> > connected_components;
  cc_extractor.Extract(&connected_components);

  // Find the largest connected component.
  int max_cc_size = 0;
  ViewId largest_cc_root_id = kInvalidViewId;
  for (const auto& connected_component : connected_components) {
    if (connected_component.second.size() > max_cc_size) {
      max_cc_size = connected_component.second.size();
      largest_cc_root_id = connected_component.first;
    }
  }

  // Remove all view pairs containing a view to remove (i.e. the ones that are
  // not in the largest connectedcomponent).
  const int num_view_pairs_before_filtering = view_graph->NumEdges();
  for (const auto& connected_component : connected_components) {
    if (connected_component.first == largest_cc_root_id) {
      continue;
    }

    // NOTE: The connected component will contain the root id as well, so we do
    // not explicity have to remove connected_component.first since it will
    // exist in connected_components.second
    for (const ViewId view_id2 : connected_component.second) {
      view_graph->RemoveView(view_id2);
      removed_views.insert(view_id2);
    }
  }

  const int num_removed_view_pairs =
      num_view_pairs_before_filtering - view_graph->NumEdges();
  LOG_IF(INFO, num_removed_view_pairs > 0)
      << num_removed_view_pairs
      << " view pairs were disconnected from the largest connected component "
         "of the view graph and were removed.";
  return removed_views;
}

}  // namespace theia
