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

#ifndef THEIA_MATCHING_FEATURE_CORRESPONDENCE_H_
#define THEIA_MATCHING_FEATURE_CORRESPONDENCE_H_

#include <cereal/access.hpp>
#include <Eigen/Core>

#include "theia/alignment/alignment.h"
#include "theia/sfm/feature.h"

namespace theia {

// The feature location of two correspondences. These can be pixel coordinates
// or normalized coordinates.
struct FeatureCorrespondence {
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  Feature feature1;
  Feature feature2;
  bool operator==(const FeatureCorrespondence& other) const {
    return (feature1 == other.feature1 && feature2 == other.feature1);
  }

 private:
  // Templated method for disk I/O with cereal. This method tells cereal which
  // data members should be used when reading/writing to/from disk.
  friend class cereal::access;
  template <class Archive>
  void serialize(Archive& ar) {  // NOLINT
    ar(feature1, feature2);
  }
};

}  // namespace theia

#endif  // THEIA_MATCHING_FEATURE_CORRESPONDENCE_H_
