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

#ifndef THEIA_SFM_EXTRACT_MAXIMALLY_PARALLEL_RIGID_SUBGRAPH_H_
#define THEIA_SFM_EXTRACT_MAXIMALLY_PARALLEL_RIGID_SUBGRAPH_H_

#include <Eigen/Core>
#include <unordered_map>

#include "theia/sfm/types.h"

namespace theia {

class ViewGraph;

// Extracts the maximally parallel rigid component of the subgraph based on
// known camera orientations and relative translation measurements.
//
// Given a set of relative translations, a rigid component is a subgraph that is
// well-posed such that global node positions may be recovered. For relative
// translations, this basically means that a triangle constraint exists on the
// node. Any nodes that are free to move (with scale or translation) with
// respect to the rigid component are not part of the rigid component. The goal
// of this method is to extract the largest rigid component so that we may
// obtain a well-posed graph for global position estimation.
//
// This method was proposed in "Identifying Maximal Rigid Components in
// Bearing-Based Localization" by Kennedy et al (IROS 2012), and was later
// utilized in "Robust Camera Location Estimation by Convex Programming" by
// Ozyesil and Singer (CVPR 2015) as a filter prior to robust global position
// estimation. Please cite these papers if using this method.
void ExtractMaximallyParallelRigidSubgraph(
    const std::unordered_map<ViewId, Eigen::Vector3d>& orientations,
    ViewGraph* view_graph);

}  // namespace theia

#endif  // THEIA_SFM_EXTRACT_MAXIMALLY_PARALLEL_RIGID_SUBGRAPH_H_
