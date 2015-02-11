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

#include <Eigen/Core>
#include <vector>
#include "gtest/gtest.h"

#include "theia/sfm/track.h"
#include "theia/sfm/types.h"
#include "theia/util/map_util.h"

namespace theia {

TEST(Track, Default) {
  Track track;
  EXPECT_EQ(track.NumViews(), 0);
  EXPECT_EQ(track.Point(), Eigen::Vector4d::Zero());
}

TEST(Track, Estimated) {
  Track track;
  EXPECT_TRUE(!track.IsEstimated());
  track.SetEstimated(true);
  EXPECT_TRUE(track.IsEstimated());
  track.SetEstimated(false);
  EXPECT_TRUE(!track.IsEstimated());
}

TEST(Track, Views) {
  Track track;
  const std::vector<ViewId> view_ids = {0, 1, 2};
  // Test that no views exist.
  EXPECT_EQ(track.NumViews(), 0);

  // Add views.
  for (int i = 0; i < view_ids.size(); i++) {
    track.AddView(view_ids[i]);
    EXPECT_EQ(track.NumViews(), i + 1);
  }

  for (int i = 0; i < view_ids.size(); i++) {
    track.RemoveView(view_ids[i]);
    EXPECT_EQ(track.NumViews(), 3 - i - 1);
  }
}

TEST(Track, ViewIds) {
  Track track;
  const std::vector<ViewId> view_ids = {0, 1, 2};

  // Add views.
  for (int i = 0; i < view_ids.size(); i++) {
    track.AddView(view_ids[i]);
  }

  // Make sure the track ids are equivalent.
  std::unordered_set<ViewId> temp_view_ids = track.ViewIds();
  EXPECT_EQ(temp_view_ids.size(), 3);
  for (int i = 0; i < view_ids.size(); i++) {
    EXPECT_TRUE(ContainsKey(temp_view_ids, view_ids[i]));
  }
}

}  // namespace theia
