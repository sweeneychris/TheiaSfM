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

#ifndef THEIA_SFM_VIEW_GRAPH_VIEW_GRAPH_H_
#define THEIA_SFM_VIEW_GRAPH_VIEW_GRAPH_H_

#include <unordered_map>
#include <unordered_set>

#include "theia/util/hash.h"
#include "theia/sfm/twoview_info.h"
#include "theia/sfm/types.h"

namespace theia {

// An undirected graph containing views in an SfM reconstruction. The graph is
// efficienctly created by only holding view ids at the vertices and
// TwoViewInfos for edge values.
//
// The implementation for an undirected view graph was insprired from:
// https://github.com/thunghan/coursera-courses/blob/master/UC%20Santa%20Cruz%20-%20C%2B%2B%20For%20C%20Programmers/Homework%202/UndirectedGraph.hpp
class ViewGraph {
 public:
  ViewGraph() {}

  // Number of views in the graph.
  int NumViews() const;

  // Number of undirected edges in the graph.
  int NumEdges() const;

  bool HasView(const ViewId view_id) const;

  // Returns a set of the ViewIds contained in the view graph.
  std::unordered_set<ViewId> ViewIds() const;

  // Removes the view from the view graph and removes all edges connected to the
  // view. Returns true on success and false if the view did not exist in the
  // view graph.
  bool RemoveView(const ViewId view_id);

  // Adds an edge between the two views with the edge value of
  // two_view_info. New vertices are added to the graph if they did not already
  // exist. If an edge already existed between the two views then the edge value
  // is updated.
  void AddEdge(const ViewId view_id_1, const ViewId view_id_2,
               const TwoViewInfo& two_view_info);

  // Removes the edge from the view graph. Returns true if the edge is removed
  // and false if the edge did not exist.
  bool RemoveEdge(const ViewId view_id_1, const ViewId view_id_2);

  // Returns the neighbor view ids for a given view, or nullptr if the view does
  // not exist.
  const std::unordered_set<ViewId>* GetNeighborIdsForView(
      const ViewId view_id) const;

  // Returns the edge value or NULL if it does not exist.
  const TwoViewInfo* GetEdge(const ViewId view_id_1,
                             const ViewId view_id_2) const;

  // Returns a map of all edges. Each edge is found exactly once in the map and
  // is indexed by the ViewIdPair (view id 1, view id 2) such that view id 1 <
  // view id 2.
  const std::unordered_map<ViewIdPair, TwoViewInfo>& GetAllEdges() const;

 private:
  // The underlying adjacency map. ViewIds are the vertices which are mapped to
  // a collection of its neighbors and the edges themselves are stored
  // separately.
  std::unordered_map<ViewId, std::unordered_set<ViewId> > vertices_;
  std::unordered_map<ViewIdPair, TwoViewInfo> edges_;
};

}  // namespace theia

#endif  // THEIA_SFM_VIEW_GRAPH_VIEW_GRAPH_H_
