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

namespace theia {

class Reconstruction;

struct BundleAdjustmentOptions {
  // For larger problems (> 1000 cameras) it is recommended to use the
  // ITERATIVE_SCHUR solver.
  ceres::LinearSolverType linear_solver_type = ceres::SPARSE_SCHUR;
  ceres::PreconditionerType preconditioner_type = ceres::SCHUR_JACOBI;
  ceres::VisibilityClusteringType visibility_clustering_type =
      ceres::SINGLE_LINKAGE;

  // If true, ceres will log verbosely.
  bool verbose = false;

  // If set to true, the camera intrinsics are held constant. This is useful if
  // the calibration is precisely known ahead of time.
  bool constant_camera_intrinsics = false;

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

}  // namespace theia

#endif  // THEIA_SFM_BUNDLE_ADJUSTMENT_BUNDLE_ADJUSTMENT_H_
