// Copyright (C) 2013 The Regents of the University of California (Regents).
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

#ifndef THEIA_SOLVERS_MLE_QUALITY_MEASUREMENT_H_
#define THEIA_SOLVERS_MLE_QUALITY_MEASUREMENT_H_

#include <glog/logging.h>
#include <algorithm>
#include <limits>
#include <memory>
#include <vector>

#include "theia/solvers/quality_measurement.h"

namespace theia {
// Define the quality metric according to Guided MLE from "MLESAC: A new robust
// estimator with application to estimating image geometry" by Torr.
class MLEQualityMeasurement : public QualityMeasurement {
 public:
  explicit MLEQualityMeasurement(const double error_thresh)
      : QualityMeasurement(error_thresh) {}

  ~MLEQualityMeasurement() {}

  // Given the residuals, assess a quality metric for the data. Returns the
  // quality assessment and outputs a vector of bools indicating the inliers.
  double ComputeCost(const std::vector<double>& residuals,
                     std::vector<int>* inliers) override {
    inliers->reserve(residuals.size());
    double mle_score = 0.0;
    for (int i = 0; i < residuals.size(); i++) {
      if (residuals[i] < error_thresh_) {
        mle_score += residuals[i];
        inliers->emplace_back(i);
      } else {
        mle_score += error_thresh_;
      }
    }
    return mle_score;
  }
};

}  // namespace theia

#endif  // THEIA_SOLVERS_MLE_QUALITY_MEASUREMENT_H_
