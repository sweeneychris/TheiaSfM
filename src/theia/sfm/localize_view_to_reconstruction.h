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

#ifndef THEIA_SFM_LOCALIZE_VIEW_TO_RECONSTRUCTION_H_
#define THEIA_SFM_LOCALIZE_VIEW_TO_RECONSTRUCTION_H_

#include "theia/sfm/bundle_adjustment/bundle_adjustment.h"
#include "theia/sfm/types.h"
#include "theia/solvers/sample_consensus_estimator.h"

namespace theia {

class Reconstruction;

// The reprojection_error_threshold_pixels is the threshold (measured in pixels)
// that determines inliers and outliers during RANSAC. This value will override
// the error thresh set in the RansacParameters.
struct LocalizeViewToReconstructionOptions {
  // The reprojection error threshold that determines whether a 2D-3D
  // correspondence is an inlier during localization.
  //
  // NOTE: This threshold is with respect to an image that is 1024 pixels
  // wide. If the image dimensions are larger or smaller than this value then
  // the threshold will be appropriately scaled.
  double reprojection_error_threshold_pixels = 4.0;

  // If true, a simplified pose solver will be used to estimate the camera
  // position given the known orientation. If that solver is not successful,
  // then standard P3P is used.
  bool assume_known_orientation = false;

  // The RANSAC parameters used for robust estimation in the localization
  // algorithms.
  RansacParameters ransac_params;

  // The view will be bundle adjusted (while all tracks are held constant) if
  // this is set to true.
  bool bundle_adjust_view = true;
  BundleAdjustmentOptions ba_options;

  // The minimum number of inliers found from RANSAC in order to be considered
  // successful localization.
  int min_num_inliers = 30;
};

// Localizes a view to the reconstruction using 2D-3D correspondences to
// estimate the absolute camera pose. If the focal length is known then the P3P
// algorithm is used, otherwise P4Pf is used to additionally recover the focal
// length.
bool LocalizeViewToReconstruction(
    const ViewId view_to_localize,
    const LocalizeViewToReconstructionOptions options,
    Reconstruction* reconstruction,
    RansacSummary* summary);

}  // namespace theia

#endif  // THEIA_SFM_LOCALIZE_VIEW_TO_RECONSTRUCTION_H_
