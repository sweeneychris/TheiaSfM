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

#ifndef THEIA_MATH_DISTRIBUTION_H_
#define THEIA_MATH_DISTRIBUTION_H_

#include <glog/logging.h>
#include <stdio.h>
#include <cmath>

#include "theia/math/util.h"

namespace theia {
// Abstract class for probability disributions.
class Distribution {
 public:
  Distribution() {}
  virtual ~Distribution() {}

  // Evaluate the disribution at x. NOTE: If eval results in a change of private
  // variables of a derived class, you should implement the changes in a
  // different method (i.e. an Update method)
  virtual double Eval(const double x) const = 0;
};

// Normal Gaussian Distribution.
class NormalDistribution : public Distribution {
 public:
  NormalDistribution(const double mean, const double sigma) : mean_(mean) {
    CHECK_GT(sigma, 0)
        << "Sigma must be greater than zero in a normal distribution";
    alpha_ = 1.0 / (sigma * sqrt(2.0 * M_PI));
    beta_ = -1.0 / (2.0 * sigma * sigma);
  }

  ~NormalDistribution() {}

  double Eval(const double x) const {
    const double normalized_x = x - mean_;
    return alpha_ * exp(beta_ * normalized_x * normalized_x);
  }

 private:
  // Normal factor.
  double alpha_;
  // Normal factor.
  double beta_;
  // The mean of the distribution.
  double mean_;
};

// Uniform distribution between left and right. Probability is uniform when x is
// within this span, and 0 when x is outside of the span.
class UniformDistribution : public Distribution {
 public:
  UniformDistribution(const double left, const double right)
      : left_(left), right_(right) {
    CHECK_LT(left, right) << "Left bound must be less than the right bound for "
                          << "uniform distributions.";
    CHECK_NE(right, left) << "Left bound is equal to right bound! Uniform "
                          << "distribution must have a nonzero range.";
    inverse_span_ = (left == right) ? 1.0 : 1.0 / (right - left);
  }

  // Destructor
  ~UniformDistribution() {}

  double Eval(const double x) const {
    return (left_ <= x && x <= right_) ? inverse_span_ : 0;
  }

 protected:
  const double left_;
  const double right_;
  double inverse_span_;
};

}  // namespace theia

#endif  // THEIA_MATH_DISTRIBUTION_H_
