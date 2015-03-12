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

#include "theia/sfm/find_common_views_by_name.h"

#include <glog/logging.h>
#include <algorithm>
#include <vector>

#include "gtest/gtest.h"
#include "theia/sfm/reconstruction.h"
#include "theia/sfm/types.h"
#include "theia/util/stringprintf.h"

namespace theia {

void TestFindCommonViewsByName(const int num_views1,
                               const int num_views2,
                               const int num_common_views) {
  Reconstruction reconstruction1, reconstruction2;
  // Add the views with the same name first.
  for (int i = 0; i < num_common_views; i++) {
    const std::string name = StringPrintf("%d", i);
    CHECK_NE(reconstruction1.AddView(name), kInvalidViewId);
    CHECK_NE(reconstruction2.AddView(name), kInvalidViewId);
  }

  // Add unique views for reconstruction1.
  for (int i = 0; i < num_views1 - num_common_views; i++) {
    const std::string name = StringPrintf("%d", reconstruction1.NumViews());
    CHECK_NE(reconstruction1.AddView(name), kInvalidViewId);
  }

  // Add unique views for reconstruction2.
  for (int i = 0; i < num_views2 - num_common_views; i++) {
    const std::string name = StringPrintf(
        "%d", reconstruction1.NumViews() + reconstruction2.NumViews());
    CHECK_NE(reconstruction2.AddView(name), kInvalidViewId);
  }

  // Validate the result.
  std::vector<std::string> common_view_names =
      FindCommonViewsByName(reconstruction1, reconstruction2);
  std::sort(common_view_names.begin(), common_view_names.end());
  EXPECT_EQ(common_view_names.size(), num_common_views);
  for (int i = 0; i < num_common_views; i++) {
    const std::string name = StringPrintf("%d", i);
    EXPECT_EQ(common_view_names[i], name);
  }
}

TEST(FindCommonViewsByName, AllViewsInCommon) {
  TestFindCommonViewsByName(10, 10, 10);
}

TEST(FindCommonViewsByName, SomeViewsInCommon) {
  TestFindCommonViewsByName(10, 10, 6);
}

TEST(FindCommonViewsByName, NoViewsInCommon) {
  TestFindCommonViewsByName(10, 10, 0);
}

}  // namespace theia
