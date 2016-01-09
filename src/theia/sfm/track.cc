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

#include "theia/sfm/track.h"

#include <Eigen/Core>
#include <unordered_set>

#include "theia/util/map_util.h"

namespace theia {

using Eigen::Vector4d;

Track::Track()
    : is_estimated_(false) {
  point_.setZero();
  color_.setZero();
}

int Track::NumViews() const {
  return view_ids_.size();
}

void Track::SetEstimated(const bool is_estimated) {
  is_estimated_ = is_estimated;
}

bool Track::IsEstimated() const {
  return is_estimated_;
}

const Eigen::Vector4d& Track::Point() const {
  return point_;
}

Vector4d* Track::MutablePoint() {
  return &point_;
}

const Eigen::Matrix<uint8_t, 3, 1>& Track::Color() const {
  return color_;
}

Eigen::Matrix<uint8_t, 3, 1>* Track::MutableColor() {
  return &color_;
}

void Track::AddView(const ViewId view_id) {
  view_ids_.insert(view_id);
}

bool Track::RemoveView(const ViewId view_id) {
  return view_ids_.erase(view_id) > 0;
}

const std::unordered_set<ViewId>& Track::ViewIds() const {
  return view_ids_;
}

}  // namespace theia
