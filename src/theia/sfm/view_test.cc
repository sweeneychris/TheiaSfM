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

#include <algorithm>
#include <string>
#include <vector>

#include "gtest/gtest.h"

#include "theia/sfm/feature.h"
#include "theia/sfm/types.h"
#include "theia/sfm/view.h"

namespace theia {

TEST(View, Name) {
  const std::string kName = "0";
  View view(kName);
  EXPECT_EQ(view.Name(), kName);
}

TEST(View, CopyConstructor) {
  static const double kView1FocalLength = 500.0;
  static const double kView2FocalLength = 1500.0;
  View view1("view_0");
  view1.MutableCamera()->SetFocalLength(kView1FocalLength);
  const std::vector<TrackId> track_ids = {0, 1, 2};
  const std::vector<Feature> features = {Feature(0, 0),
                                         Feature(1, 1),
                                         Feature(2, 2)};
  // Add features.
  for (int i = 0; i < track_ids.size(); i++) {
    view1.AddFeature(track_ids[i], features[i]);
  }

  // Use the copy constructor and test that some internal variables are equal.
  View view2(view1);
  EXPECT_EQ(view1.Camera().FocalLength(), view2.Camera().FocalLength());
  EXPECT_EQ(view1.NumFeatures(), view2.NumFeatures());

  // Alter the view2 camera parameters and ensure this does not affect view1's
  // camera parameters.
  view2.MutableCamera()->SetFocalLength(kView2FocalLength);
  EXPECT_NE(view1.Camera().FocalLength(), view2.Camera().FocalLength());
}

TEST(View, Estimated) {
  View view;
  EXPECT_TRUE(!view.IsEstimated());
  view.SetEstimated(true);
  EXPECT_TRUE(view.IsEstimated());
  view.SetEstimated(false);
  EXPECT_TRUE(!view.IsEstimated());
}

TEST(View, Features) {
  View view;
  const std::vector<TrackId> track_ids = {0, 1, 2};
  const std::vector<Feature> features = {Feature(0, 0),
                                         Feature(1, 1),
                                         Feature(2, 2)};
  // Test that no features exist.
  EXPECT_EQ(view.NumFeatures(), 0);

  // Add features.
  for (int i = 0; i < track_ids.size(); i++) {
    view.AddFeature(track_ids[i], features[i]);
    const Feature* feature = view.GetFeature(track_ids[i]);
    EXPECT_NE(feature, nullptr);
    EXPECT_EQ(*feature, features[i]);
    EXPECT_EQ(view.NumFeatures(), i + 1);
  }

  for (int i = 0; i < track_ids.size(); i++) {
    view.RemoveFeature(track_ids[i]);
    const Feature* feature = view.GetFeature(track_ids[i]);
    EXPECT_EQ(feature, nullptr);
    EXPECT_EQ(view.NumFeatures(), 3 - i - 1);
  }
}


TEST(View, TrackIds) {
  View view;
  const std::vector<TrackId> track_ids = {0, 1, 2};
  const std::vector<Feature> features = {Feature(0, 0),
                                         Feature(1, 1),
                                         Feature(2, 2)};
  // Add features.
  for (int i = 0; i < track_ids.size(); i++) {
    view.AddFeature(track_ids[i], features[i]);
  }

  // Make sure the track ids are equivalent.
  std::vector<TrackId> temp_track_ids = view.TrackIds();
  std::sort(temp_track_ids.begin(), temp_track_ids.end());
  for (int i = 0; i < track_ids.size(); i++) {
    EXPECT_EQ(track_ids[i], temp_track_ids[i]);
  }
}

}  // namespace theia
