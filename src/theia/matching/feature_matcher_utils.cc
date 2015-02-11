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

#include "theia/matching/feature_matcher_utils.h"

#include <glog/logging.h>
#include <unordered_map>
#include <vector>

#include "theia/matching/indexed_feature_match.h"
#include "theia/util/map_util.h"

namespace theia {

// Modifies forward matches so that it removes all matches that are not
// contained in the backwards matches.
void IntersectMatches(const std::vector<IndexedFeatureMatch>& backwards_matches,
                      std::vector<IndexedFeatureMatch>* forward_matches) {
  std::unordered_map<int, int> index_map;
  index_map.reserve(backwards_matches.size());
  // Add all feature2 -> feature1 matches to the map.
  for (const IndexedFeatureMatch& feature_match : backwards_matches) {
    InsertOrDie(&index_map,
                feature_match.feature1_ind,
                feature_match.feature2_ind);
  }

  // Search the map for feature1 -> feature2 matches that are also present in
  // the feature2 -> feature1 matches.
  auto match_iterator = forward_matches->begin();
  while (match_iterator != forward_matches->end()) {
    if (match_iterator->feature1_ind !=
        FindWithDefault(index_map, match_iterator->feature2_ind, -1)) {
      match_iterator = forward_matches->erase(match_iterator);
      continue;
    }

    ++match_iterator;
  }
}

}  // namespace theia
