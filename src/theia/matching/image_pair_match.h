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

#ifndef THEIA_MATCHING_IMAGE_PAIR_MATCH_H_
#define THEIA_MATCHING_IMAGE_PAIR_MATCH_H_

#include <cereal/access.hpp>
#include <cereal/types/vector.hpp>
#include <vector>

#include "theia/alignment/alignment.h"
#include "theia/matching/feature_correspondence.h"
#include "theia/sfm/twoview_info.h"

namespace theia {

struct ImagePairMatch {
 public:
  std::string image1;
  std::string image2;

  // If the matches are verified matches then the two view info contains the
  // relative pose information between the images.
  TwoViewInfo twoview_info;

  // Feature locations in pixel coordinates. If the match is a verified match
  // then this only contains inlier correspondences.
  std::vector<FeatureCorrespondence> correspondences;

 private:
  // Templated method for disk I/O with cereal. This method tells cereal which
  // data members should be used when reading/writing to/from disk.
  friend class cereal::access;
  template <class Archive>
  void serialize(Archive& ar) {  // NOLINT
    ar(image1, image2, twoview_info, correspondences);
  }
};

struct ImagePairMatchDeprecated {
 public:
  // Indices of the matches image pair with respect to the input vectors.
  int image1_index;
  int image2_index;

  // If the matches are verified matches then the two view info contains the
  // relative pose information between the images.
  TwoViewInfo twoview_info;

  // Feature locations in pixel coordinates. If the match is a verified match
  // then this only contains inlier correspondences.
  std::vector<FeatureCorrespondence> correspondences;

 private:
  // Templated method for disk I/O with cereal. This method tells cereal which
  // data members should be used when reading/writing to/from disk.
  friend class cereal::access;
  template <class Archive>
  void serialize(Archive& ar) {  // NOLINT
    ar(image1_index, image2_index, twoview_info, correspondences);
  }
};

}  // namespace theia

#endif  // THEIA_MATCHING_IMAGE_PAIR_MATCH_H_
