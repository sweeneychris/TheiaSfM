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

#ifndef THEIA_SFM_BUNDLE_ADJUSTMENT_BUNDLE_ADJUSTMENT_H_
#define THEIA_SFM_BUNDLE_ADJUSTMENT_BUNDLE_ADJUSTMENT_H_

#include <ceres/ceres.h>
#include <unordered_set>

#include "theia/sfm/bundle_adjustment/create_loss_function.h"
#include "theia/sfm/types.h"
#include "theia/util/enable_enum_bitmask_operators.h"

namespace theia {

class Reconstruction;

// The camera intrinsics parameters are defined by:
//   - Focal length
//   - Aspect ratio
//   - Skew
//   - Principal points (x and y)
//   - Radial distortion (2-parameter model)
// It is often known for instance that skew is 0 and aspect ratio is 1, and so
// we do not always desire to optimize all camera intrinsics. In many cases, the
// focal length is the only parameter we care to optimize.
//
// Users can specify which intrinsics to optimize by using a bitmask. For
// instance FOCAL_LENGTH|PRINCIPAL_POINTS will optimize the focal length and
// principal points. The options NONE and ALL are also given for convenience.
enum class OptimizeIntrinsicsType {
  NONE = 0x00,
  FOCAL_LENGTH = 0x01,
  ASPECT_RATIO = 0x02,
  SKEW = 0x04,
  PRINCIPAL_POINTS = 0x08,
  RADIAL_DISTORTION = 0x10,
  ALL =
    FOCAL_LENGTH | ASPECT_RATIO | SKEW | PRINCIPAL_POINTS | RADIAL_DISTORTION,
};
ENABLE_ENUM_BITMASK_OPERATORS(OptimizeIntrinsicsType)

struct BundleAdjustmentOptions {
  // The type of loss function used for BA. By default, we use a standard L2
  // loss function, but robust cost functions could be used.
  LossFunctionType loss_function_type = LossFunctionType::TRIVIAL;
  double robust_loss_width = 10.0;

  // For larger problems (> 1000 cameras) it is recommended to use the
  // ITERATIVE_SCHUR solver.
  ceres::LinearSolverType linear_solver_type = ceres::SPARSE_SCHUR;
  ceres::PreconditionerType preconditioner_type = ceres::SCHUR_JACOBI;
  ceres::VisibilityClusteringType visibility_clustering_type =
      ceres::CANONICAL_VIEWS;

  // If true, ceres will log verbosely.
  bool verbose = false;

  // Indicates which intrinsics should be optimized as part of bundle
  // adjustment. By default, we do not optimize skew and aspect ratio since
  // these are almost universally constant.
  OptimizeIntrinsicsType intrinsics_to_optimize =
      OptimizeIntrinsicsType::FOCAL_LENGTH |
      OptimizeIntrinsicsType::PRINCIPAL_POINTS |
      OptimizeIntrinsicsType::RADIAL_DISTORTION;

  int num_threads = 1;
  int max_num_iterations = 500;

  // Max BA time is 1 hour.
  double max_solver_time_in_seconds = 3600.0;

  // Inner iterations can improve the quality according to the Ceres email list.
  bool use_inner_iterations = true;

  // These variables may be useful to change if the optimization is converging
  // to a bad result.
  double function_tolerance = 1e-6;
  double gradient_tolerance = 1e-10;
  double parameter_tolerance = 1e-8;
  double max_trust_region_radius = 1e12;
};

// Some important metrics for analyzing bundle adjustment results.
struct BundleAdjustmentSummary {
  // This only indicates whether the optimization was successfully run and makes
  // no guarantees on the quality or convergence.
  bool success = false;
  double initial_cost = 0.0;
  double final_cost = 0.0;
  double setup_time_in_seconds = 0.0;
  double solve_time_in_seconds = 0.0;
};

// Bundle adjust all views and tracks in the reconstruction.
BundleAdjustmentSummary BundleAdjustReconstruction(
    const BundleAdjustmentOptions& options, Reconstruction* reconstruction);

// Bundle adjust the specified views and all tracks observed by those views.
BundleAdjustmentSummary BundleAdjustPartialReconstruction(
    const BundleAdjustmentOptions& options,
    const std::unordered_set<ViewId>& views_to_optimize,
    const std::unordered_set<TrackId>& tracks_to_optimize,
    Reconstruction* reconstruction);

// Bundle adjust a single view.
BundleAdjustmentSummary BundleAdjustView(const BundleAdjustmentOptions& options,
                                         const ViewId view_id,
                                         Reconstruction* reconstruction);

// Bundle adjust a single track.
BundleAdjustmentSummary BundleAdjustTrack(const BundleAdjustmentOptions& options,
                                          const TrackId track_id,
                                          Reconstruction* reconstruction);

}  // namespace theia

#endif  // THEIA_SFM_BUNDLE_ADJUSTMENT_BUNDLE_ADJUSTMENT_H_
