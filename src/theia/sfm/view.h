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

#ifndef THEIA_SFM_VIEW_H_
#define THEIA_SFM_VIEW_H_

#include <cereal/access.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/unordered_map.hpp>
#include <string>
#include <unordered_map>
#include <vector>

#include "theia/sfm/camera/camera.h"
#include "theia/sfm/camera_intrinsics_prior.h"
#include "theia/sfm/feature.h"
#include "theia/sfm/types.h"

namespace theia {

// A View contains high level information about an image that has been
// captured. This includes the name, EXIF metadata, and track information that
// is found through feature matching.
class View {
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  View();
  explicit View(const std::string& name);

  ~View() {}

  const std::string& Name() const;

  void SetEstimated(bool is_estimated);
  bool IsEstimated() const;

  const class Camera& Camera() const;
  class Camera* MutableCamera();

  const struct CameraIntrinsicsPrior& CameraIntrinsicsPrior() const;
  struct CameraIntrinsicsPrior* MutableCameraIntrinsicsPrior();

  int NumFeatures() const;

  std::vector<TrackId> TrackIds() const;

  const Feature* GetFeature(const TrackId track_id) const;

  void AddFeature(const TrackId track_id, const Feature& feature);

  bool RemoveFeature(const TrackId track_id);

 private:
  // Templated method for disk I/O with cereal. This method tells cereal which
  // data members should be used when reading/writing to/from disk.
  friend class cereal::access;
  template <class Archive>
  void serialize(Archive& ar) {  // NOLINT
    ar(name_, is_estimated_, camera_, camera_intrinsics_prior_, features_);
  }

  std::string name_;
  bool is_estimated_;
  class Camera camera_;
  struct CameraIntrinsicsPrior camera_intrinsics_prior_;
  std::unordered_map<TrackId, Feature> features_;
};

}  // namespace theia

#endif  // THEIA_SFM_VIEW_H_
