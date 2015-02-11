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

#include "theia/math/sturm_chain.h"

#include <Eigen/Core>
#include <glog/logging.h>
#include <algorithm>
#include <cmath>
#include <vector>

#include "theia/math/polynomial.h"

namespace theia {

using Eigen::VectorXd;

namespace {

int CountSignChanges(const std::vector<bool>& vals) {
  int num_sign_changes = 0;
  for (int i = 1; i < vals.size(); i++) {
    if (vals[i] != vals[i - 1]) {
      num_sign_changes++;
    }
  }
  return num_sign_changes;
}

}  // namespace

SturmChain::SturmChain(const Eigen::VectorXd& polynomial) {
  sturm_chain_.reserve(polynomial.size());
  sturm_chain_.push_back(polynomial);
  sturm_chain_.push_back(DifferentiatePolynomial(polynomial));
  while (sturm_chain_.back().size() > 1) {
    const VectorXd& poly1 = sturm_chain_[sturm_chain_.size() - 2];
    const VectorXd& poly2 = sturm_chain_[sturm_chain_.size() - 1];
    VectorXd remainder1 = poly1.tail(poly1.size() - 1);
    remainder1.head(remainder1.size() - 1) -=
        poly1(0) / poly2(0) * poly2.tail(poly2.size() - 1);

    const VectorXd remainder2 =
        remainder1.tail(remainder1.size() - 1) -
        remainder1(0) / poly2(0) * poly2.tail(poly2.size() - 1);
    sturm_chain_.emplace_back(-remainder2);
  }
}

int SturmChain::NumSignChanges(const double x) const {
  std::vector<bool> sturm_chain_sign(sturm_chain_.size());
  for (int i = 0; i < sturm_chain_.size(); i++) {
    sturm_chain_sign[i] = EvaluatePolynomial(sturm_chain_[i], x) > 0;
  }
  return CountSignChanges(sturm_chain_sign);
}

void SturmChain::ComputeRootBounds(double* lower_bound,
                                   double* upper_bound) {
  const auto& poly = sturm_chain_[0];
  const int n = poly.size();
  const double sqrt_part =
      sqrt(poly[1] * poly[1] - (2.0 * n * poly[0] * poly[2]) / (n - 1.0));
  const double first_part = -poly[1] / (n * poly[0]);
  const double second_part = (n - 1.0) / (n * poly[0]);
  *lower_bound = first_part - second_part * sqrt_part;
  *upper_bound = first_part + second_part * sqrt_part;
  if (*lower_bound > *upper_bound) {
    std::swap(*lower_bound, *upper_bound);
  }
}

int SturmChain::NumSignChangesAtInfinity() const {
  std::vector<bool> pos_infinity_sign(sturm_chain_.size());
  for (int i = 0; i < sturm_chain_.size(); i++) {
    pos_infinity_sign[i] = sturm_chain_[i][0] > 0;
  }
  return CountSignChanges(pos_infinity_sign);
}

int SturmChain::NumSignChangesAtNegativeInfinity() const {
  std::vector<bool> neg_infinity_sign(sturm_chain_.size());
  bool odd_degree = sturm_chain_.size() % 2 == 0;
  for (int i = 0; i < sturm_chain_.size(); i++) {
    // Negative infinty evaluates to negative if the leading coeff is negative
    // and the degree is even or the leading coeff is positive and the
    // degree is odd.
    neg_infinity_sign[i] = (sturm_chain_[i][0] < 0) != odd_degree;
    odd_degree = !odd_degree;
  }

  return CountSignChanges(neg_infinity_sign);
}

}  // namespace theia
