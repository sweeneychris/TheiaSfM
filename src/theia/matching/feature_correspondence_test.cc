// Copyright (C) 2016 The Regents of the University of California (Regents).
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

#include <glog/logging.h>
#include "gtest/gtest.h"

#include "theia/matching/feature_correspondence.h"

namespace theia {

TEST(FeatureCorrespondence, Equality) {
  FeatureCorrespondence correspondence;
  correspondence.feature1 = Eigen::Vector2d(1.0, 2.0);
  correspondence.feature2 = Eigen::Vector2d(3.0, 4.0);

  FeatureCorrespondence correspondence2 = correspondence;
  EXPECT_TRUE(correspondence == correspondence2);
}

TEST(FeatureCorrespondence, Inequality) {
  FeatureCorrespondence correspondence;
  correspondence.feature1 = Eigen::Vector2d(1.0, 2.0);
  correspondence.feature2 = Eigen::Vector2d(3.0, 4.0);

  FeatureCorrespondence correspondence2;
  correspondence2.feature1 = Eigen::Vector2d(3.0, 4.0);
  correspondence2.feature2 = Eigen::Vector2d(1.0, 2.0);
  EXPECT_FALSE(correspondence == correspondence2);
}

}  // namespace theia
