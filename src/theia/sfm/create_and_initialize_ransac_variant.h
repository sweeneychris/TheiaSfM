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

#ifndef THEIA_SFM_CREATE_AND_INITIALIZE_RANSAC_VARIANT_H_
#define THEIA_SFM_CREATE_AND_INITIALIZE_RANSAC_VARIANT_H_

#include <glog/logging.h>

#include "theia/solvers/lmed.h"
#include "theia/solvers/prosac.h"
#include "theia/solvers/sample_consensus_estimator.h"
#include "theia/solvers/ransac.h"

namespace theia {

// NOTE: Prosac requires correspondences to be sorted by the descriptor
// distances with the best match first. See theia/solvers for more information
// on the various types.
enum class RansacType {
  RANSAC = 0,
  PROSAC = 1,
  LMED = 2
};

// Factory method to create a ransac variant based on the specified options. The
// variante is then initialized and fails if initialization is not successful.
template <class Estimator>
std::unique_ptr<SampleConsensusEstimator<Estimator> >
CreateAndInitializeRansacVariant(
    const RansacType& ransac_type,
    const RansacParameters& ransac_options, const Estimator& estimator) {
  std::unique_ptr<SampleConsensusEstimator<Estimator> > ransac_variant;
  switch (ransac_type) {
    case RansacType::RANSAC:
      ransac_variant.reset(new Ransac<Estimator>(ransac_options, estimator));
      break;
    case RansacType::PROSAC:
      ransac_variant.reset(new Prosac<Estimator>(ransac_options, estimator));
      break;
    case RansacType::LMED:
      ransac_variant.reset(new LMed<Estimator>(ransac_options, estimator));
      break;
    default:
      ransac_variant.reset(new Ransac<Estimator>(ransac_options, estimator));
      break;
  }

  CHECK(ransac_variant->Initialize()) << "Could not initialize ransac "
                                         "estimator for estimating two view "
                                         "reconstructions";
  return ransac_variant;
}

}  // namespace theia

#endif  // THEIA_SFM_CREATE_AND_INITIALIZE_RANSAC_VARIANT_H_
