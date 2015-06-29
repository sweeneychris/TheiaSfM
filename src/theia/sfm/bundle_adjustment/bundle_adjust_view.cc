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

#include "theia/sfm/bundle_adjustment/bundle_adjust_view.h"

#include <ceres/ceres.h>
#include <vector>

#include "theia/sfm/camera/reprojection_error.h"
#include "theia/sfm/reconstruction.h"
#include "theia/sfm/track.h"
#include "theia/sfm/types.h"
#include "theia/sfm/view.h"

namespace theia {
namespace {
void AddCameraParametersToProblem(const bool constant_camera_intrinsics,
                                  double* camera_parameters,
                                  ceres::Problem* problem) {
  if (constant_camera_intrinsics) {
    std::vector<int> constant_intrinsics;
    for (int i = 0; i < Camera::kIntrinsicsSize; i++) {
      constant_intrinsics.push_back(Camera::kExtrinsicsSize + i);
    }
    ceres::SubsetParameterization* subset_parameterization =
        new ceres::SubsetParameterization(Camera::kParameterSize,
                                          constant_intrinsics);
    problem->AddParameterBlock(camera_parameters,
                               Camera::kParameterSize,
                               subset_parameterization);
  } else {
    problem->AddParameterBlock(camera_parameters, Camera::kParameterSize);
  }

  // Set bounds for certain camera parameters to make sure they are reasonable.
  problem->SetParameterLowerBound(
      camera_parameters, Camera::kExtrinsicsSize + Camera::FOCAL_LENGTH, 0.0);
  problem->SetParameterLowerBound(
      camera_parameters, Camera::kExtrinsicsSize + Camera::ASPECT_RATIO, 0.0);
}

}  // namespace

// Bundle adjust a track. This will keep all views constant and only modify the
// track position.
bool BundleAdjustView(const ViewId view_id,
                      const bool optimize_camera_intrinsics,
                      Reconstruction* reconstruction) {
  CHECK_NOTNULL(reconstruction);
  View* view = CHECK_NOTNULL(reconstruction->MutableView(view_id));
  if (!view->IsEstimated()) {
    return false;
  }

  ceres::Problem problem;

  // Set solver options.
  ceres::Solver::Options solver_options;
  solver_options.logging_type = ceres::SILENT;
  solver_options.linear_solver_type = ceres::DENSE_QR;

  // Add the camera parameters to the problem. Optionally set the intrinsics to
  // constant.
  Camera* camera = view->MutableCamera();
  AddCameraParametersToProblem(optimize_camera_intrinsics,
                               camera->mutable_parameters(),
                               &problem);

  const auto& track_ids = view->TrackIds();
  for (const TrackId track_id : track_ids) {
    Track* track = CHECK_NOTNULL(reconstruction->MutableTrack(track_id));
    // Only optimize over tracks that are estimated.
    if (!track->IsEstimated()) {
      continue;
    }

    const Feature* feature = CHECK_NOTNULL(view->GetFeature(track_id));
    problem.AddResidualBlock(
        ReprojectionError::Create(*feature),
        NULL,
        camera->mutable_parameters(),
        track->MutablePoint()->data());
    problem.SetParameterBlockConstant(track->MutablePoint()->data());
  }

  // Solve the problem.
  ceres::Solver::Summary solver_summary;
  ceres::Solve(solver_options, &problem, &solver_summary);

  return solver_summary.IsSolutionUsable();
}

}  // namespace theia
