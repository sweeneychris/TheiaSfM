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
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVERCAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#ifndef STATX_DISTRIBUTIONS_RAYLEIGH_H_
#define STATX_DISTRIBUTIONS_RAYLEIGH_H_

#include <cmath>
#include <vector>

namespace statx {
namespace distributions {

// Calculates the Rayleigh density at x with sigma as its parameter
inline double raylpdf(const double x,
                      const double sigma) {
  const double sigma_sqrd = sigma*sigma;
  return x * exp(-0.5 * x*x / sigma_sqrd) / sigma_sqrd;
}

// Calculates the Rayleigh cummulative probability at x with sigma as its
// parameter.
inline double raylcdf(const double x,
                      const double sigma) {
  const double sigma_sqrd = sigma*sigma;
  return 1.0 - exp(-0.5 * x*x / sigma_sqrd);
}

// Computes the maximum likelihood (ML) estimate of the Rayleigh distribution.
inline double raylfit(const std::vector<double>& samples) {
  double sigma = 0.0;
  for (double z : samples) sigma += z*z;
  sigma = 0.5 * sigma / samples.size();
  return sqrt(sigma);
}
}  // distributions
}  // statx
#endif  // STATX_DISTRIBUTIONS_RAYLEIGH_H_
