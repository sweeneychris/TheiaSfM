// Copyright (C) 2014  Victor Fragoso <vfragoso@cs.ucsb.edu>
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
//     * Neither the name of the University of California, Santa Barbara nor the
//       names of its contributors may be used to endorse or promote products
//       derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL VICTOR FRAGOSO BE LIABLE FOR ANY DIRECT,
// INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#ifndef STATX_DISTRIBUTIONS_GAMMA_H_
#define STATX_DISTRIBUTIONS_GAMMA_H_

#include <cmath>
#include <vector>
#include <limits>

namespace statx {
namespace distributions {
/// Evaluates the lower incomplete gamma
/// Parameters a & x >= 0
inline double lower_inc_gamma(const double a,
                              const double x) {
  if (a <= 0.0 || x <= 0.0) return std::numeric_limits<double>::infinity();
  double acc = 0.0;
  double term = 1.0 / a;
  int n = 1;
  while (term != 0) {
    acc += term;
    term *= x / (a + n++);
  }
  return pow(x, a)*exp(-x)*acc;
}

// Evaluates the gamma density at x with k and theta parameters.
// Parameter constraints:
// x > 0
// k > 0
// theta > 0
inline double gammapdf(const double x,
                       const double k,
                       const double theta) {
  if ( x <= 0 || k <= 0 || theta <= 0) return 0.0;
  double term1 = 1.0 / (tgamma(k) * pow(theta, k));
  double term2 = exp(-x / theta);
  return term1 * pow(x, k - 1) * term2;
}

// Evaluates the gamma distribution function (CDF) at x with k and theta
// parameters.
// Parameter constraints:
// x > 0
// k > 0
// theta > 0
inline double gammacdf(const double x,
                       const double k,
                       const double theta) {
  if (x <= 0.0 || k <= 0.0 || theta <= 0.0) return 0.0;
  double term1 = 1.0 / tgamma(k);
  double term2 = lower_inc_gamma(k, x / theta);
  return term1*term2;
}

// Computes the Maximum Likelihood estimate of the Gamma distribution
// k: shape parameter
// theta: scale parameter
bool gammafit(const std::vector<double>& data,
              double* k,  // shape
              double* theta /* scale*/);
}  // distributions
}  // statx
#endif  // STATX_DISTRIBUTIONS_GAMMA_H_
