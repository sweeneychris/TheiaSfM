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

#ifndef THEIA_SFM_FILTER_VIEW_PAIRS_FROM_ORIENTATION_H_
#define THEIA_SFM_FILTER_VIEW_PAIRS_FROM_ORIENTATION_H_

#include <Eigen/Core>
#include <unordered_map>

#include "theia/util/hash.h"
#include "theia/sfm/twoview_info.h"
#include "theia/sfm/types.h"

namespace theia {

// Filters view pairs based on the orientation estimates. If the relative
// rotation obtained from the two view match (i.e. TwoViewInfo.rotation_2)
// differs from the relative rotation formed by the orientation estimates by
// more than the threshold then the view pair is deemed "bad" and is removed.
// The view id pairs asssumes the view with the lower id be the first entry in
// the ViewIdPair.
//
// Concretely, we only keep R_{i,j} if: ||R_{i,j} - R_j * R_i^t|| < tau, where
// || A - B || is the angular distance between A and B.
//
// NOTE: This function will remove any view pairs that contain a view that does
// not have an entry in orientations and will output a comment to LOG(WARNING).
void FilterViewPairsFromOrientation(
    const std::unordered_map<ViewId, Eigen::Vector3d>& orientations,
    const double max_relative_rotation_difference_degrees,
    std::unordered_map<ViewIdPair, TwoViewInfo>* view_pairs);

}  // namespace theia

#endif  // THEIA_SFM_FILTER_VIEW_PAIRS_FROM_ORIENTATION_H_
