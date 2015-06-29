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

#include "theia/sfm/filter_view_pairs_from_relative_translation.h"

#include <ceres/rotation.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <memory>
#include <mutex>  // NOLINT
#include <unordered_map>

#include "theia/math/util.h"
#include "theia/util/hash.h"
#include "theia/util/map_util.h"
#include "theia/util/random.h"
#include "theia/util/threadpool.h"
#include "theia/sfm/twoview_info.h"
#include "theia/sfm/types.h"
#include "theia/sfm/view_graph/view_graph.h"

namespace theia {
using Eigen::Vector3d;

namespace {

// Helper struct to maintain the graph for the translation projection problem.
struct MFASNode {
  std::unordered_map<ViewId, double> incoming_nodes;
  std::unordered_map<ViewId, double> outgoing_nodes;
  double incoming_weight = 0;
  double outgoing_weight = 0;
};

// Rotate the translation direction based on the known orientation such that the
// translation is in the global reference frame.
std::unordered_map<ViewIdPair, Vector3d>
      RotateRelativeTranslationsToGlobalFrame(
          const std::unordered_map<ViewId, Vector3d>& orientations,
          const std::unordered_map<ViewIdPair, TwoViewInfo>& view_pairs) {
  std::unordered_map<ViewIdPair, Vector3d> rotated_translations;
  rotated_translations.reserve(orientations.size());

  for (const auto& view_pair : view_pairs) {
    const Vector3d view_to_world_rotation =
        -1.0 * FindOrDie(orientations, view_pair.first.first);
    Vector3d rotated_translation;
    ceres::AngleAxisRotatePoint(view_to_world_rotation.data(),
                                view_pair.second.position_2.data(),
                                rotated_translation.data());
    rotated_translations.emplace(view_pair.first, rotated_translation);
  }
  return rotated_translations;
}

// Find the next view to add to the order. We attempt to choose a source (i.e.,
// a node with no incoming edges) or choose a node based on a heuristic such
// that it has the most source-like properties.
ViewId FindNextViewInOrder(
    const std::unordered_map<ViewId, MFASNode>& degrees_for_view) {
  ViewId best_choice = kInvalidViewId;
  double best_score = 0;
  for (const auto& view : degrees_for_view) {
    // If the view is a source view, return it.
    if (view.second.incoming_nodes.size() == 0) {
      return view.first;
    }

    // Otherwise, keep track of the max score seen so far.
    const double score = (view.second.outgoing_weight + 1.0) /
                         (view.second.incoming_weight + 1.0);
    if (score > best_score) {
      best_choice = view.first;
      best_score = score;
    }
  }

  return best_choice;
}

// Based on the 1D translation projections, compute an ordering of the
// translations.
std::unordered_map<ViewId, int> OrderTranslationsFromProjections(
    const std::unordered_map<ViewIdPair, double>&
        translation_direction_projections) {
  // Compute the degrees of all vertices as the sum of weights coming in or out.
  std::unordered_map<ViewId, MFASNode> degrees_for_view;
  for (const auto& translation_projection : translation_direction_projections) {
    const ViewIdPair view_id_pair =
        (translation_projection.second > 0)
            ? translation_projection.first
            : ViewIdPair(translation_projection.first.second,
                         translation_projection.first.first);

    // Update the MFAS entry.
    const double weight = std::abs(translation_projection.second);
    degrees_for_view[view_id_pair.second].incoming_weight += weight;
    degrees_for_view[view_id_pair.first].outgoing_weight += weight;
    degrees_for_view[view_id_pair.second].incoming_nodes.emplace(
        view_id_pair.first, weight);
    degrees_for_view[view_id_pair.first].outgoing_nodes.emplace(
        view_id_pair.second, weight);
  }

  // Compute the ordering.
  const int num_views = degrees_for_view.size();
  std::unordered_map<ViewId, int> translation_ordering;
  for (int i = 0; i < num_views; i++) {
    // Find the next view to add.
    const ViewId next_view_in_order = FindNextViewInOrder(degrees_for_view);
    translation_ordering[next_view_in_order] = i;

    // Update the MFAS graph and remove the next view from the degrees_for_view.
    const auto& next_view_info =
        FindOrDie(degrees_for_view, next_view_in_order);
    for (auto& neighbor_info : next_view_info.incoming_nodes) {
      degrees_for_view[neighbor_info.first].outgoing_weight -=
          neighbor_info.second;
      degrees_for_view[neighbor_info.first].outgoing_nodes.erase(
          next_view_in_order);
    }
    for (auto& neighbor_info : next_view_info.outgoing_nodes) {
      degrees_for_view[neighbor_info.first].incoming_weight -=
          neighbor_info.second;
      degrees_for_view[neighbor_info.first].incoming_nodes.erase(
          next_view_in_order);
    }
    degrees_for_view.erase(next_view_in_order);
  }

  return translation_ordering;
}

// Projects all the of the translation onto the given axis.
std::unordered_map<ViewIdPair, double> ProjectTranslationsOntoAxis(
    const Vector3d& axis,
    const std::unordered_map<ViewIdPair, Vector3d>& relative_translations) {
  std::unordered_map<ViewIdPair, double> projection_weights;
  projection_weights.reserve(relative_translations.size());

  for (const auto& relative_translation : relative_translations) {
    const double projection_weight = relative_translation.second.dot(axis);
    projection_weights.emplace(relative_translation.first, projection_weight);
  }
  return projection_weights;
}

// This chooses a random axis based on the given relative translations.
void ComputeMeanVariance(
    const std::unordered_map<ViewIdPair, Vector3d>& relative_translations,
    Vector3d* mean,
    Vector3d* variance) {
  mean->setZero();
  variance->setZero();
  for (const auto& translation : relative_translations) {
    *mean += translation.second;
  }
  *mean /= static_cast<double>(relative_translations.size());

  for (const auto& translation : relative_translations) {
    *variance += (translation.second - *mean).cwiseAbs2();
  }
  *variance /= static_cast<double>(relative_translations.size() - 1);
}

// Performs a single iterations of the translation filtering. This method is
// thread-safe.
void TranslationFilteringIteration(
    const std::unordered_map<ViewIdPair, Vector3d>& relative_translations,
    const Vector3d& direction_mean,
    const Vector3d& direction_variance,
    std::mutex* mutex,
    std::unordered_map<ViewIdPair, double>* bad_edge_weight) {
  // Get a random vector to project all relative translations on to.
  const Vector3d random_axis =
      Vector3d(RandGaussian(direction_mean[0], direction_variance[0]),
               RandGaussian(direction_mean[1], direction_variance[1]),
               RandGaussian(direction_mean[2], direction_variance[2]))
          .normalized();

  // Project all vectors.
  const std::unordered_map<ViewIdPair, double>&
      translation_direction_projections =
          ProjectTranslationsOntoAxis(random_axis, relative_translations);

  // Compute ordering.
  const std::unordered_map<ViewId, int>& translation_ordering =
      OrderTranslationsFromProjections(translation_direction_projections);

  // Compute bad edge weights.
  for (auto& edge : *bad_edge_weight) {
    const int ordering_diff =
        FindOrDie(translation_ordering, edge.first.second) -
        FindOrDie(translation_ordering, edge.first.first);
    const double& projection_weight_of_edge =
        FindOrDieNoPrint(translation_direction_projections, edge.first);

    VLOG(3) << "Edge (" << edge.first.first << ", " << edge.first.second
            << ") has ordering diff of " << ordering_diff
            << " and a projection of " << projection_weight_of_edge << " from "
            << FindOrDieNoPrint(relative_translations, edge.first).transpose();
    // If the ordering is inconsistent, add the absolute value of the bad weight
    // to the aggregate bad weight.
    if ((ordering_diff < 0 && projection_weight_of_edge > 0) ||
        (ordering_diff > 0 && projection_weight_of_edge < 0)) {
      mutex->lock();
      edge.second += std::abs(projection_weight_of_edge);
      mutex->unlock();
    }
  }
}

}  // namespace

void FilterViewPairsFromRelativeTranslation(
    const FilterViewPairsFromRelativeTranslationOptions& options,
    const std::unordered_map<ViewId, Vector3d>& orientations,
    ViewGraph* view_graph) {
  const auto& view_pairs = view_graph->GetAllEdges();

  // Weights of edges that have been accumulated throughout the iterations. A
  // higher weight means the edge is more likely to be bad.
  std::unordered_map<ViewIdPair, double> bad_edge_weight;
  for (const auto& view_pair : view_pairs) {
    bad_edge_weight[view_pair.first] = 0.0;
  }

  // Compute the adjusted translations so that they are oriented in the global
  // frame.
  const std::unordered_map<ViewIdPair, Vector3d>& rotated_translations =
      RotateRelativeTranslationsToGlobalFrame(orientations, view_pairs);

  Vector3d translation_mean, translation_variance;
  ComputeMeanVariance(rotated_translations,
                      &translation_mean,
                      &translation_variance);

  std::unique_ptr<ThreadPool> pool(new ThreadPool(options.num_threads));
  std::mutex mutex;
  for (int i = 0; i < options.num_iterations; i++) {
    pool->Add(TranslationFilteringIteration,
              rotated_translations,
              translation_mean,
              translation_variance,
              &mutex,
              &bad_edge_weight);
  }
  pool.reset(nullptr);

  // Remove all the bad edges.
  const double max_aggregated_projection_tolerance =
      options.translation_projection_tolerance * options.num_iterations;
  int num_view_pairs_removed = 0;
  for (const auto& view_pair : bad_edge_weight) {
    VLOG(3) << "View pair (" << view_pair.first.first << ", "
            << view_pair.first.second << ") projection = " << view_pair.second;
    if (view_pair.second > max_aggregated_projection_tolerance) {
      view_graph->RemoveEdge(view_pair.first.first, view_pair.first.second);
      ++num_view_pairs_removed;
    }
  }

  VLOG(1) << "Removed " << num_view_pairs_removed
          << " view pairs by relative translation filtering.";
}

}  // namespace theia
