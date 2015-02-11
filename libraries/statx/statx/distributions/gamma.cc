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

#include "statx/distributions/gamma.h"
#include <cmath>
#include <asa121.hpp>
#include <asa103.hpp>

namespace statx {
namespace distributions {
using std::vector;
using asa103::digamma;
using asa121::trigamma;

// Implements the Maximum Likelihood estimation algorithm.
// The MLE problem for the Gamma distribution is a concave one.
// The implementation maximizes the Likelihood following an
// iterative scheme. We follow the algorithm described in Wikipedia.
bool gammafit(const vector<double>& data,
              double* k,  // shape
              double* theta /* scale*/) {
  if (data.empty() || !k || !theta) return false;

  // Compute the initial value of s
  double acc1 = 0.0;
  double acc2 = 0.0;
  for (double z : data) {
    acc1 += z;
    acc2 += log(z);
  }

  int fault;
  const double epsilon = 1e-6;
  const int max_iters = 100;
  const int n = data.size();
  const double mean = acc1/n;
  const double s = log(mean) - acc2/n;
  const double s_term = s - 3;
  *k = (3 - s + sqrt(s_term * s_term + 24 * s)) / (12 * s);

  // Updating loop
  for (int i = 0; i < max_iters; i++) {
    const double k_inv = 1.0 / *k;
    const double digamma_val = digamma(*k, &fault);
    if (fault) return false;
    const double trigamma_val = trigamma(*k, &fault);
    if (fault) return false;
    const double delta_k =
        (log(*k) - digamma_val - s) / (k_inv - trigamma_val);
    *k -= delta_k;
    if (delta_k < epsilon) break;
  }

  *theta = mean / (*k);
  return true;
}
}  // distributions
}  // statx
