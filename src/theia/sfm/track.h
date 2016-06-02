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

#ifndef THEIA_SFM_TRACK_H_
#define THEIA_SFM_TRACK_H_

#include <cereal/access.hpp>
#include <cereal/cereal.hpp>
#include <cereal/types/unordered_set.hpp>
#include <Eigen/Core>
#include <stdint.h>
#include <unordered_set>

#include "theia/io/eigen_serializable.h"
#include "theia/sfm/types.h"

namespace theia {

// A track contains information about a 3D point and the views that observe the
// point. This is based off of LibMV's Structure class:
// https://github.com/libmv/libmv/blob/master/src/libmv/multiview/structure.h
class Track {
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  Track();
  ~Track() {}

  int NumViews() const;

  void SetEstimated(const bool is_estimated);
  bool IsEstimated() const;

  const Eigen::Vector4d& Point() const;
  Eigen::Vector4d* MutablePoint();

  const Eigen::Matrix<uint8_t, 3, 1>& Color() const;
  Eigen::Matrix<uint8_t, 3, 1>* MutableColor();

  void AddView(const ViewId view_id);
  bool RemoveView(const ViewId view_id);

  const std::unordered_set<ViewId>& ViewIds() const;

 private:
  // Templated method for disk I/O with cereal. This method tells cereal which
  // data members should be used when reading/writing to/from disk.
  friend class cereal::access;
  template <class Archive>
  void serialize(Archive& ar, const std::uint32_t version) {  // NOLINT
    ar(is_estimated_, view_ids_, point_, color_);
  }

  bool is_estimated_;
  std::unordered_set<ViewId> view_ids_;
  Eigen::Vector4d point_;
  Eigen::Matrix<uint8_t, 3, 1> color_;
};

}  // namespace theia

CEREAL_CLASS_VERSION(theia::Track, 0);

#endif  // THEIA_SFM_TRACK_H_
