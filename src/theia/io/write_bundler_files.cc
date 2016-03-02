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

#include "theia/io/write_bundler_files.h"

#include <Eigen/Core>
#include <glog/logging.h>
#include <fstream>  // NOLINT
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "theia/sfm/camera/camera.h"
#include "theia/sfm/camera_intrinsics_prior.h"
#include "theia/sfm/reconstruction_estimator_utils.h"
#include "theia/sfm/reconstruction.h"
#include "theia/sfm/types.h"
#include "theia/sfm/track.h"
#include "theia/sfm/view.h"
#include "theia/util/map_util.h"

namespace theia {
namespace {

bool WriteBundleFile(const Reconstruction& reconstruction,
                     const std::string& bundle_file,
                     const std::string& lists_file) {
  // Output file stream for bundle file.
  std::ofstream ofs_bundle(bundle_file);
  if (!ofs_bundle.is_open()) {
    LOG(ERROR) << "Cannot open the file: " << bundle_file << " for writing.";
    return false;
  }

  // Output file stream for list file.
  std::ofstream ofs_lists(lists_file);
  if (!ofs_lists.is_open()) {
    LOG(ERROR) << "Cannot open the file: " << lists_file << " for writing.";
    return false;
  }

  const Eigen::Matrix3d theia_to_bundler =
      Eigen::Vector3d(1.0, -1.0, -1.0).asDiagonal();

  ofs_bundle << "# Bundle file v0.3" << std::endl;
  ofs_bundle << reconstruction.NumViews() << " " << reconstruction.NumTracks()
      << std::endl;

  const Eigen::IOFormat unaligned(Eigen::FullPrecision, Eigen::DontAlignCols);
  // Output all cameras first.
  std::unordered_map<ViewId, int> view_id_to_index;
  const auto& view_ids = reconstruction.ViewIds();
  for (int i = 0; i < view_ids.size(); i++) {
    view_id_to_index[view_ids[i]] = i;
    const View* view = reconstruction.View(view_ids[i]);

    // Output the information to the list file.
    ofs_lists << view->Name();
    const auto& prior = view->CameraIntrinsicsPrior();
    if (prior.focal_length.is_set) {
      ofs_lists << " 0 " << prior.focal_length.value;
    }
    ofs_lists << std::endl;

    // Output the camera information to the bundle file.
    const Camera& camera = view->Camera();
    ofs_bundle << camera.FocalLength() << " " << camera.RadialDistortion1()
               << " " << camera.RadialDistortion2() << std::endl;

    const Eigen::Matrix3d rotation =
        theia_to_bundler * camera.GetOrientationAsRotationMatrix();
    ofs_bundle << rotation.format(unaligned) << std::endl;

    const Eigen::Vector3d translation =
        theia_to_bundler *
        (-camera.GetOrientationAsRotationMatrix() * camera.GetPosition());
    ofs_bundle << translation.transpose().format(unaligned) << std::endl;
  }

  // Output all points
  const auto& track_ids = reconstruction.TrackIds();
  for (const TrackId track_id : track_ids) {
    const Track* track = reconstruction.Track(track_id);
    const Eigen::Vector3d position = track->Point().hnormalized();
    ofs_bundle << position.transpose().format(unaligned) << std::endl;

    ofs_bundle << track->Color().cast<double>()[0] << " "
               << track->Color().cast<double>()[1] << " "
               << track->Color().cast<double>()[2] << std::endl;
    const auto& views_in_track = track->ViewIds();
    ofs_bundle << views_in_track.size();
    for (const ViewId view_id : views_in_track) {
      const int index = FindOrDie(view_id_to_index, view_id);
      const View* view = reconstruction.View(view_id);
      const Feature* feature = view->GetFeature(track_id);
      // Bundler has pixel coordinates with the origin at the center of the
      // image and positive x to the right, positive y is up.
      const Eigen::Vector2d adjusted_feature(
          feature->x() - view->Camera().PrincipalPointX(),
          -(feature->y() - view->Camera().PrincipalPointY()));

      // Note we give the keypoint index as 0 because we do not store SIFT
      // keyfiles.
      ofs_bundle << " " << index << " 0 "
          << adjusted_feature.transpose().format(unaligned);
    }
    ofs_bundle << std::endl;
  }
  return true;
}

}  // namespace

bool WriteBundlerFiles(const Reconstruction& reconstruction,
                       const std::string& lists_file,
                       const std::string& bundle_file) {
  Reconstruction estimated_reconstruction;
  CreateEstimatedSubreconstruction(reconstruction, &estimated_reconstruction);

  return WriteBundleFile(estimated_reconstruction, bundle_file, lists_file);
}

}  // namespace theia
