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

#include "theia/sfm/colorize_reconstruction.h"

#include <Eigen/Core>

#include <mutex>  // NOLINT
#include <string>
#include <memory>

#include "theia/image/image.h"
#include "theia/sfm/reconstruction.h"
#include "theia/sfm/track.h"
#include "theia/sfm/types.h"
#include "theia/sfm/view.h"
#include "theia/util/filesystem.h"
#include "theia/util/map_util.h"
#include "theia/util/threadpool.h"

namespace theia {
namespace {

void ExtractColorsFromImage(
    const std::string& image_file,
    const View& view,
    std::unordered_map<TrackId, Eigen::Vector3f>* colors,
    std::mutex* mutex_lock) {
  LOG(INFO) << "Extracting color for features in image: " << image_file;
  const FloatImage image(image_file);

  const auto& track_ids = view.TrackIds();
  if (image.Channels() == 3) {
    for (const TrackId track_id : track_ids) {
      const Feature feature = *view.GetFeature(track_id);
      const int x = static_cast<int>(feature.x());
      const int y = static_cast<int>(feature.y());
      std::lock_guard<std::mutex> lock(*mutex_lock);
      (*colors)[track_id] += 255.0 * image.GetXY(x, y);
    }
  } else if (image.Channels() == 1) {
    for (const TrackId track_id : track_ids) {
      const Feature feature = *view.GetFeature(track_id);
      const int x = static_cast<int>(feature.x());
      const int y = static_cast<int>(feature.y());
      std::lock_guard<std::mutex> lock(*mutex_lock);
      (*colors)[track_id] += 255.0 * image.GetXY(x, y);
    }
  } else {
    LOG(FATAL) << "The image file at: " << image_file
               << " is not an RGB or a grayscale image so the color cannot be "
                  "extracted.";
  }
}

}  // namespace

void ColorizeReconstruction(const std::string& image_directory,
                            const int num_threads,
                            Reconstruction* reconstruction) {
  CHECK(DirectoryExists(image_directory))
      << "The image directory " << image_directory << " does not exist.";
  CHECK_GT(num_threads, 0);
  CHECK_NOTNULL(reconstruction);

  ThreadPool pool(num_threads);
  std::mutex mutex_lock;

  // Initialize the colors to be (0, 0, 0).
  const auto& track_ids = reconstruction->TrackIds();
  std::unordered_map<TrackId, Eigen::Vector3f> colors;
  for (const TrackId track_id : track_ids) {
    colors[track_id].setZero();
  }

  // For each image, find the color of each feature and add the value to the
  // colors map.
  const auto& view_ids = reconstruction->ViewIds();
  for (const ViewId view_id : view_ids) {
    const View* view = reconstruction->View(view_id);
    const std::string image_filepath = image_directory + view->Name();
    CHECK(FileExists(image_filepath)) << "The image file: " << image_filepath
                                      << " does not exist!";
    pool.Add(ExtractColorsFromImage,
             image_filepath,
             *view,
             &colors,
             &mutex_lock);
  }
  pool.WaitForTasksToFinish();

  // The colors map now contains a sum of all colors, so to get the mean we must
  // divide by the number of observations in each track.
  for (const TrackId track_id : track_ids) {
    Eigen::Vector3f color = FindOrDie(colors, track_id);

    Track* track = reconstruction->MutableTrack(track_id);
    color /= static_cast<float>(track->NumViews());
    *track->MutableColor() = color.cast<uint8_t>();
  }
}

}  // namespace theia
