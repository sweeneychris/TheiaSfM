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

#ifndef THEIA_MATCHING_FEATURE_MATCHER_OPTIONS_H_
#define THEIA_MATCHING_FEATURE_MATCHER_OPTIONS_H_

#include <string>

#include "theia/sfm/verify_two_view_matches.h"

namespace theia {

// Options for matching image collections.
struct FeatureMatcherOptions {
  // Number of threads to use in parallel for matching.
  int num_threads = 1;

  // Matching may be performed in core (i.e. all in memory) or out-of-core. For
  // the latter, features are written and read to/from disk as needed (utilizing
  // an LRU cache). The out-of-core strategy is more scalable since the memory
  // footprint is limited. Set this value to false to perform all-in-memory
  // matching.
  bool match_out_of_core = false;

  // Keypoints and descriptors are stored to disk as they are added to the
  // FeatureMatcher. Features will be stored in this directory, which must be a
  // valid writeable directory.
  std::string keypoints_and_descriptors_output_dir = "";

  // We store the descriptors of up to cache_capacity images in the cache at a
  // given time. The higher the cache capacity, the more memory is required to
  // perform image-to-image matching.
  int cache_capacity = 128;

  // Only symmetric matches are kept.
  bool keep_only_symmetric_matches = true;

  // Only keep the matches that pass the lowes ratio test such that the distance
  // of the best best match is less than lowes_ratio of the distance of the
  // second nearest neighbor match.
  bool use_lowes_ratio = true;
  float lowes_ratio = 0.8;

  // After performing feature matching with descriptors typically the 2-view
  // geometry is estimated using RANSAC (from the matched descriptors) and only
  // the features that support the estimated geometry are "verified" as
  // plausible and are kept. If set to true, geometric verification will be
  // performed to obtain higher quality matches.
  bool perform_geometric_verification = true;

  // The parameter settings for geometric verification.
  VerifyTwoViewMatchesOptions geometric_verification_options;

  // Only images that contain more feature matches than this number will be
  // returned.
  int min_num_feature_matches = 30;
};

}  // namespace theia

#endif  // THEIA_MATCHING_FEATURE_MATCHER_OPTIONS_H_
