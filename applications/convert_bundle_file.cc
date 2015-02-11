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

#include <Eigen/Core>
#include <glog/logging.h>
#include <gflags/gflags.h>
#include <theia/theia.h>

#include <algorithm>
#include <memory>
#include <string>


DEFINE_string(lists_file, "", "Input bundle lists file.");
DEFINE_string(bundle_file, "", "Input bundle file.");
DEFINE_string(output_reconstruction_file, "",
              "Output reconstruction file in binary format.");

int main(int argc, char* argv[]) {
  google::InitGoogleLogging(argv[0]);
  google::ParseCommandLineFlags(&argc, &argv, true);

  // Load the SIFT descriptors into the cameras.
  std::unique_ptr<theia::Reconstruction> reconstruction(
      new theia::Reconstruction());
  CHECK(theia::ReadBundlerFiles(FLAGS_lists_file,
                                FLAGS_bundle_file,
                                reconstruction.get()))
      << "Could not read Bundler files.";

  // Check that the reprojection errors are sane.
  int num_projections_behind_camera = 0;
  std::vector<double> reprojection_errors;
  for (const theia::TrackId track_id : reconstruction->TrackIds()) {
    const theia::Track* track = CHECK_NOTNULL(reconstruction->Track(track_id));
    for (const theia::ViewId view_id : track->ViewIds()) {
      const theia::Feature* feature =
          reconstruction->View(view_id)->GetFeature(track_id);

      // Reproject the observations.
      Eigen::Vector2d projection;
      if (reconstruction->View(view_id)
              ->Camera().ProjectPoint(track->Point(), &projection) < 0) {
        ++num_projections_behind_camera;
      }

      // Compute reprojection error.
      const double reprojection_error = (*feature - projection).norm();
      reprojection_errors.emplace_back(reprojection_error);
    }
  }

  std::sort(reprojection_errors.begin(), reprojection_errors.end());
  const double mean_reprojection_error =
      std::accumulate(reprojection_errors.begin(),
                      reprojection_errors.end(),
                      0.0) / static_cast<double>(reprojection_errors.size());
  const double median_reprojection_error =
      reprojection_errors[reprojection_errors.size() / 2];

  LOG(INFO) << "\nNum views: " << reconstruction->NumViews()
            << "\nNum 3D points: " << reconstruction->NumTracks()
            << "\nNum observations: " << reprojection_errors.size()
            << "\nNum reprojections behind camera: "
            << num_projections_behind_camera
            << "\nMean reprojection error = " << mean_reprojection_error
            << "\nMedian reprojection_error = " << median_reprojection_error;

  CHECK(WriteReconstruction(*reconstruction, FLAGS_output_reconstruction_file))
      << "Could not write out reconstruction file.";
  return 0;
}
