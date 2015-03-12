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

#ifndef THEIA_SFM_RECONSTRUCTION_ESTIMATOR_H_
#define THEIA_SFM_RECONSTRUCTION_ESTIMATOR_H_

#include <unordered_set>

#include "theia/sfm/types.h"
#include "theia/sfm/reconstruction_estimator_options.h"

namespace theia {

class Reconstruction;
struct ReconstructionEstimatorOptions;
class ViewGraph;

struct ReconstructionEstimatorSummary {
  bool success = false;
  std::unordered_set<ViewId> estimated_views;
  std::unordered_set<TrackId> estimated_tracks;

  // All times are given in seconds.
  double initial_view_graph_filtering_time = 0.0;
  double camera_intrinsics_calibration_time = 0.0;
  double rotation_estimation_time = 0.0;
  double rotation_filtering_time = 0.0;
  double relative_translation_optimization_time = 0.0;
  double relative_translation_filtering_time = 0.0;
  double position_estimation_time = 0.0;
  double triangulation_time = 0.0;
  double bundle_adjustment_time = 0.0;
};

// A reconstruction estimator should build a reconstruction from a view graph
// and an unestimated reconstruction with the corresponding views of the view
// graph. The camera poses and 3D point positions are estimated with this
// class. The focus of the Theia library is to create global methods for SfM
// such that all camera pose are estimated simultaneously then points are
// estimated afterwards.
class ReconstructionEstimator {
 public:
  virtual ~ReconstructionEstimator() {}

  // Estimates the camera poses for a reconstruction given the view graph
  // describing the multi-view correspondences.
  virtual ReconstructionEstimatorSummary Estimate(
      ViewGraph* view_graph, Reconstruction* reconstruction) = 0;

  static ReconstructionEstimator* Create(
      const ReconstructionEstimatorOptions& options);
};

}  // namespace theia

#endif  // THEIA_SFM_RECONSTRUCTION_ESTIMATOR_H_
