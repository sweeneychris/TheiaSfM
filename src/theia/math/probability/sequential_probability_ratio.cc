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

#include "theia/math/probability/sequential_probability_ratio.h"

#include <cmath>

namespace theia {

double CalculateSPRTDecisionThreshold(double sigma, double epsilon,
                                      double time_compute_model_ratio,
                                      int num_models_verified) {
  // Eq. 2 in Matas et. al.
  double c = (1.0 - sigma) * log((1.0 - sigma) / (1.0 - epsilon)) +
             sigma * log(sigma / epsilon);

  // Eq. 6 in Matas et. al.
  double a_0 = time_compute_model_ratio * c /
                   static_cast<double>(num_models_verified) + 1.0;
  double decision_threshold = a_0;
  double kConvergence = 1e-4;
  // Matas et. al. says the decision threshold typically converges in 4
  // iterations. Set the max iters to 1000 as a safeguard, but test for
  // convergence.
  for (int i = 0; i < 1000; i++) {
    double new_decision_threshold = a_0 + log(decision_threshold);
    double step_difference = fabs(new_decision_threshold - decision_threshold);
    decision_threshold = new_decision_threshold;
    // If the series has converged, break.
    if (step_difference < kConvergence) break;
  }
  return decision_threshold;
}

bool SequentialProbabilityRatioTest(const std::vector<double>& residuals,
                                    double error_thresh, double sigma,
                                    double epsilon, double decision_threshold,
                                    int* num_tested_points,
                                    double* observed_inlier_ratio) {
  int observed_num_inliers = 0;
  double likelihood_ratio = 1.0;
  for (int i = 0; i < residuals.size(); i++) {
    // Check whether i-th data point is consistent with the model. Update the
    // likelihood ratio accordingly.
    if (residuals[i] < error_thresh) {
      likelihood_ratio *= sigma / epsilon;
      observed_num_inliers += 1;
    } else {
      likelihood_ratio *= (1.0 - sigma) / (1.0 - epsilon);
    }

    // If likehood ratio exceeds our decision threshold we can terminate early.
    if (likelihood_ratio > decision_threshold) {
      *observed_inlier_ratio = static_cast<double>(observed_num_inliers) /
                               static_cast<double>(i + 1);
      *num_tested_points = i + 1;
      return false;
    }
  }

  *observed_inlier_ratio = static_cast<double>(observed_num_inliers) /
                           static_cast<double>(residuals.size());
  *num_tested_points = residuals.size();
  return true;
}
}  // namespace theia
