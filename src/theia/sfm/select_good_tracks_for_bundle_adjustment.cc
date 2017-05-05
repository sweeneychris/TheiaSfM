// Copyright (C) 2017 The Regents of the University of California (Regents).
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

#include "theia/sfm/select_good_tracks_for_bundle_adjustment.h"

#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "theia/sfm/camera/camera.h"
#include "theia/sfm/reconstruction.h"
#include "theia/sfm/reconstruction_estimator_utils.h"
#include "theia/sfm/track.h"
#include "theia/sfm/types.h"
#include "theia/sfm/view.h"
#include "theia/util/hash.h"
#include "theia/util/map_util.h"

namespace theia {
namespace {
// Track statistics are the track length and mean reprojection error.
typedef std::pair<int, double> TrackStatistics;
typedef std::pair<TrackId, TrackStatistics> GridCellElement;
typedef std::unordered_map<Eigen::Vector2i, std::vector<GridCellElement> >
    ImageGrid;

// Sorts the grid cell elements by the track statistics, which will sort first
// by the (truncated) track length, then by the mean reprojection error.
bool CompareGridCellElements(
    const std::pair<TrackId, TrackStatistics>& element1,
    const std::pair<TrackId, TrackStatistics>& element2) {
  return element1.second < element2.second;
}

// Return the squared reprojection error of the track in the view.
inline double ComputeSqReprojectionError(const View& view,
                                         const Feature& feature,
                                         const Track& track) {
  Eigen::Vector2d reprojected_feature;
  view.Camera().ProjectPoint(track.Point(), &reprojected_feature);
  return (reprojected_feature - feature).squaredNorm();
}

// Compute the reprojection error and truncated track length for this specific
// track.
TrackStatistics ComputeStatisticsForTrack(
    const Reconstruction& reconstruction,
    const TrackId track_id,
    const int long_track_length_threshold) {
  // Any tracks that reach this function are guaranteed to exist and be
  // estimated, so no need to check for that here.
  const Track* track = reconstruction.Track(track_id);
  const auto& views_observing_track = track->ViewIds();

  double sq_reprojection_error_sum = 0.0;
  int num_valid_reprojections = 0;
  // Compute the sq reprojection error for each view that observes the track
  // and it it to the accumulating sum.
  for (const ViewId view_id : views_observing_track) {
    const View* view = reconstruction.View(view_id);
    if (view == nullptr || !view->IsEstimated()) {
      continue;
    }
    sq_reprojection_error_sum +=
        ComputeSqReprojectionError(*view, *view->GetFeature(track_id), *track);
    ++num_valid_reprojections;
  }

  // Compute and return the track statistics.
  const int truncated_track_length =
      std::min(num_valid_reprojections, long_track_length_threshold);
  const double mean_sq_reprojection_error =
      sq_reprojection_error_sum / static_cast<double>(num_valid_reprojections);
  return TrackStatistics(truncated_track_length, mean_sq_reprojection_error);
}

// Compute the mean reprojection error and the truncated track length of each
// track. We truncate the track length based on the observation that while
// larger track lengths provide better constraints for bundle adjustment, larger
// tracks are also more likely to contain outliers in our experience. Truncating
// the track lengths enforces that the long tracks with the lowest reprojection
// error are chosen.
void ComputeTrackStatistics(
    const Reconstruction& reconstruction,
    const std::unordered_set<ViewId>& view_ids,
    const int long_track_length_threshold,
    std::unordered_map<TrackId, TrackStatistics>* track_statistics) {
  // Iterate over all views and compute the track statistics for each track we
  // encounter.
  for (const ViewId view_id : view_ids) {
    const View* view = reconstruction.View(view_id);
    const auto& tracks_in_view = view->TrackIds();
    // Compute statistics for each track in this view that we have not already
    // computed statistics for.
    for (const TrackId track_id : tracks_in_view) {
      const Track* track = reconstruction.Track(track_id);
      // Skip this track if the statistics were already computed, or if the
      // track has not been estimated.
      if (ContainsKey(*track_statistics, track_id) || track == nullptr ||
          !track->IsEstimated()) {
        continue;
      }

      // Compute the track statistics and add it to the output map.
      const TrackStatistics& statistics_for_this_track =
          ComputeStatisticsForTrack(reconstruction,
                                    track_id,
                                    long_track_length_threshold);
      track_statistics->emplace(track_id, statistics_for_this_track);
    }
  }
}

// Select tracks from the image to ensure good spatial coverage of the image. To
// do this, we first bin the tracks into grid cells in an image grid. Then
// within each cell we find the best ranked track and add it to the list of
// tracks to optimize.
void SelectBestTracksFromEachImageGridCell(
    const Reconstruction& reconstruction,
    const View& view,
    const int grid_cell_size,
    const std::unordered_map<TrackId, TrackStatistics>& track_statistics,
    std::unordered_set<TrackId>* tracks_to_optimize) {
  static const double inv_grid_cell_size =  1.0 / grid_cell_size;

  // Hash each feature into a grid cell.
  ImageGrid image_grid;
  const auto& track_ids = view.TrackIds();
  for (const TrackId track_id : track_ids) {
    const Track* track = reconstruction.Track(track_id);
    if (track == nullptr || !track->IsEstimated()) {
      continue;
    }

    const Feature& feature = *view.GetFeature(track_id);
    const TrackStatistics& current_track_statistics =
        FindOrDieNoPrint(track_statistics, track_id);
    const Eigen::Vector2i grid_cell =
        (feature * inv_grid_cell_size).cast<int>();

    image_grid[grid_cell].emplace_back(track_id, current_track_statistics);
  }

  // Select the best feature from each grid cell and add it to the tracks to
  // optimize.
  for (auto& grid_cell : image_grid) {
    // Order the features in each cell by track length first, then mean
    // reprojection error.
    const GridCellElement& grid_cell_element =
        *std::min_element(grid_cell.second.begin(),
                          grid_cell.second.end(),
                          CompareGridCellElements);

    // Insert the track id in to the tracks to optimize.
    tracks_to_optimize->emplace(grid_cell_element.first);
  }
}

// Selects the top ranked tracks that have not already been chosen until the
// view observes the minimum number of optimized tracks.
void SelectTopRankedTracksInView(
    const Reconstruction& reconstruction,
    const std::unordered_map<TrackId, TrackStatistics>& track_statistics,
    const View& view,
    const int min_num_optimized_tracks_per_view,
    std::unordered_set<TrackId>* tracks_to_optimize) {
  int num_optimized_tracks = 0;
  int num_estimated_tracks = 0;

  const auto& tracks_in_view = view.TrackIds();
  std::vector<GridCellElement> ranked_candidate_tracks;
  for (const TrackId track_id : tracks_in_view) {
    const Track* track = reconstruction.Track(track_id);
    if (track == nullptr || !track->IsEstimated()) {
      continue;
    }
    // We only reach this point if the track is estimated.
    ++num_estimated_tracks;

    // If the track is already slated for optimization, increase the count of
    // optimized features.
    if (ContainsKey(*tracks_to_optimize, track_id)) {
      ++num_optimized_tracks;
      // If the number of optimized_tracks is greater than the minimum then we
      // can return early since we know that no more features need to added for
      // this view.
      if (num_optimized_tracks >= min_num_optimized_tracks_per_view) {
        return;
      }
    } else {
      // If the track is not already set to be optimized then add it to the list
      // of candidate tracks.
      ranked_candidate_tracks.emplace_back(
          track_id, FindOrDieNoPrint(track_statistics, track_id));
    }
  }

  // We only reach this point if the number of optimized tracks is less than the
  // minimum. If that is the case then we add the top candidate features until
  // the minimum number of features observed is met.
  if (num_optimized_tracks != num_estimated_tracks) {
    // Select how many tracks to add. If we need more tracks than are estimated
    // then we simply add all remaining features.
    const int num_optimized_tracks_needed =
        std::min(min_num_optimized_tracks_per_view - num_optimized_tracks,
                 num_estimated_tracks - num_optimized_tracks);
    std::partial_sort(
        ranked_candidate_tracks.begin(),
        ranked_candidate_tracks.begin() + num_optimized_tracks_needed,
        ranked_candidate_tracks.end());
    // Add the candidate tracks to the list of tracks to be optimized.
    for (int i = 0; i < num_optimized_tracks_needed; i++) {
      tracks_to_optimize->emplace(ranked_candidate_tracks[i].first);
    }
  }
}

}  // namespace

// The efficiency of large scale bundle adjustment can be dramatically increased
// by choosing only a subset of 3d points to optimize, as the 3d points tend to
// have increasing scene redundancy. If the points are chosen in a way that
// properly constrains the nonlinear optimization, similar results in quality
// may be observed compared to optimized all tracks.
//
// The 3d points are chosen such that they fit the following criteria:
//    a) High confidence (i.e. low reprojection error).
//    b) Long tracks are preferred.
//    c) The tracks used for optimization provide a good spatial coverage in
//       each image.
bool SelectGoodTracksForBundleAdjustment(
    const Reconstruction& reconstruction,
    const int long_track_length_threshold,
    const int image_grid_cell_size_pixels,
    const int min_num_optimized_tracks_per_view,
    std::unordered_set<TrackId>* tracks_to_optimize) {
  std::unordered_set<ViewId> view_ids;
  GetEstimatedViewsFromReconstruction(reconstruction, &view_ids);
  return SelectGoodTracksForBundleAdjustment(reconstruction,
                                             view_ids,
                                             long_track_length_threshold,
                                             image_grid_cell_size_pixels,
                                             min_num_optimized_tracks_per_view,
                                             tracks_to_optimize);
}

bool SelectGoodTracksForBundleAdjustment(
    const Reconstruction& reconstruction,
    const std::unordered_set<ViewId>& view_ids,
    const int long_track_length_threshold,
    const int image_grid_cell_size_pixels,
    const int min_num_optimized_tracks_per_view,
    std::unordered_set<TrackId>* tracks_to_optimize) {
  // Compute the track mean reprojection errors.
  std::unordered_map<TrackId, TrackStatistics> track_statistics;
  ComputeTrackStatistics(reconstruction,
                         view_ids,
                         long_track_length_threshold,
                         &track_statistics);

  // For each image, divide the image into a grid and choose the highest quality
  // tracks from each grid cell. This encourages good spatial coverage of tracks
  // within each image.
  for (const ViewId view_id : view_ids) {
    const View* view = reconstruction.View(view_id);

    // Select the best tracks from each grid cell in the image and add them to
    // the container of tracks to be optimized.
    SelectBestTracksFromEachImageGridCell(reconstruction,
                                          *view,
                                          image_grid_cell_size_pixels,
                                          track_statistics,
                                          tracks_to_optimize);
  }

  // To this point, we have only added features that have as full spatial
  // coverage as possible within each image but we have not ensured that each
  // image is constrainted by at least K features. So, we cycle through all
  // views again and add the top M tracks that have not already been added.
  for (const ViewId view_id : view_ids) {
    const View* view = reconstruction.View(view_id);

    // If this view is not constrained by enough optimized tracks, add the top
    // ranked features until there are enough tracks constraining the view.
    SelectTopRankedTracksInView(reconstruction,
                                track_statistics,
                                *view,
                                min_num_optimized_tracks_per_view,
                                tracks_to_optimize);
  }

  return true;
}

}  // namespace theia
