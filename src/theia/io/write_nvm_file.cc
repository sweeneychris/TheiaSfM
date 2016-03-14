// Copyright (C) 2016 The Regents of the University of California (Regents).
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

#include "theia/io/write_nvm_file.h"

#include <glog/logging.h>
#include <fstream>
#include <string>
#include <unordered_map>

#include "theia/sfm/reconstruction.h"
#include "theia/sfm/track.h"
#include "theia/sfm/view.h"
#include "theia/util/map_util.h"

namespace theia {

// Writes an NVM file that may then be inspected with Visual SfM or other
// software packages.
bool WriteNVMFile(const std::string& nvm_filepath,
                  const Reconstruction& reconstruction) {
  std::ofstream nvm_file;
  nvm_file.open(nvm_filepath.c_str());
  if (!nvm_file.is_open()) {
    LOG(WARNING) << "Could not open nvm file for writing: " << nvm_filepath;
    return false;
  }

  // Output the NVM header.
  nvm_file << "NVM_V3 " << std::endl << std::endl;

  // Number of cameras.
  const auto& view_ids = reconstruction.ViewIds();
  nvm_file << view_ids.size() << std::endl;
  std::unordered_map<ViewId, int> view_id_to_index;
  std::unordered_map<ViewId, std::unordered_map<TrackId, int> >
      feature_index_mapping;
  // Output each camera.
  for (const ViewId view_id : view_ids) {
    const int current_index = view_id_to_index.size();
    view_id_to_index[view_id] = current_index;

    const View& view = *reconstruction.View(view_id);
    const Camera& camera = view.Camera();
    const Eigen::Quaterniond quat(camera.GetOrientationAsRotationMatrix());
    const Eigen::Vector3d position(camera.GetPosition());
    nvm_file << view.Name() << " " << camera.FocalLength() << " " << quat.w()
             << " " << quat.x() << " " << quat.y() << " " << quat.z() << " "
             << position.x() << " " << position.y() << " " << position.z()
             << " " << camera.RadialDistortion1() << " 0" << std::endl;

    // Assign each feature in this view to a unique feature index (unique within
    // each image, not unique to the reconstruction).
    const auto& view_track_ids = reconstruction.View(view_id)->TrackIds();
    for (int i = 0; i < view_track_ids.size(); i++) {
      const TrackId track_id = view_track_ids[i];
      feature_index_mapping[view_id][track_id] = i;
    }
  }

  // Number of points.
  const auto& track_ids = reconstruction.TrackIds();
  nvm_file << track_ids.size() << std::endl;
  // Output each point.
  for (const TrackId track_id : track_ids) {
    const Track* track = reconstruction.Track(track_id);
    const Eigen::Vector3d position = track->Point().hnormalized();

    // Normalize the color.
    Eigen::Vector3i color = track->Color().cast<int>();

    nvm_file << position.x() << " " << position.y() << " " << position.z()
             << " " << color.x() << " " << color.y() << " " << color.z() << " "
             << track->NumViews() << " ";

    // Output the observations of this 3D point.
    const auto& views_observing_track = track->ViewIds();
    for (const ViewId& view_id : views_observing_track) {
      const View* view = reconstruction.View(view_id);

      // Get the feature location normalized by the principal point.
      const Camera& camera = view->Camera();
      const Feature feature =
          (*view->GetFeature(track_id)) -
          Feature(camera.PrincipalPointX(), camera.PrincipalPointY());

      const int track_index =
          FindOrDie(FindOrDie(feature_index_mapping, view_id), track_id);
      const int view_index = FindOrDie(view_id_to_index, view_id);
      nvm_file << view_index << " " << track_index << " " << feature.x() << " "
               << feature.y() << " ";
    }
    nvm_file << std::endl;
  }

  // Indicate the end of the file.
  nvm_file << "0" << std::endl;
  nvm_file.close();
  return true;
}

}  // namespace theia
