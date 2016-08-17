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

#ifndef THEIA_SFM_FILTER_VIEW_PAIRS_FROM_RELATIVE_TRANSLATION_H_
#define THEIA_SFM_FILTER_VIEW_PAIRS_FROM_RELATIVE_TRANSLATION_H_

#include <Eigen/Core>
#include <memory>
#include <unordered_map>

#include "theia/sfm/types.h"
#include "theia/util/random.h"

namespace theia {

class ViewGraph;

struct FilterViewPairsFromRelativeTranslationOptions {
  // Random number generator used to sample from normal distributions to compute
  // the axis of projection for 1dsfm filtering.
  std::shared_ptr<RandomNumberGenerator> rng;

  // Filtering the translations is embarassingly parallel (each iteration can
  // be run independently) so we can use a thread pool to speed up computation.
  int num_threads = 1;

  // The projection will be performed for the given number of iterations (we
  // recommend > 40 iterations).
  int num_iterations = 48;

  // The parameter translation_projection_tolerance determines which
  // translations are considered "bad" after analyzing their projections over
  // many iterations (it corresponds to tau in the paper).
  double translation_projection_tolerance = 0.08;
};

// Filters view pairs based on the relative translation estimations according to
// the algorithm presented in "Robust Global Translations with 1DSfM" by Kyle
// Wilson and Noah Snavely (ECCV 2014). This algorithm determines translations
// directions which are likely to be outliers by projecting translations to a
// 1-dimensional subproblem. The relative translations are repeatedly projected
// onto a (semi) random vector and are ordered to find a consistent
// embedding. Translation projections which are inconsistent with the ordering
// are likely to be outliers. This process is repeated over many iterations to
// determine translation directions likely to be outliers. Please see the paper
// for more specific details.
void FilterViewPairsFromRelativeTranslation(
    const FilterViewPairsFromRelativeTranslationOptions& options,
    const std::unordered_map<ViewId, Eigen::Vector3d>& orientations,
    ViewGraph* view_graph);

}  // namespace theia

#endif  // THEIA_SFM_FILTER_VIEW_PAIRS_FROM_RELATIVE_TRANSLATION_H_
