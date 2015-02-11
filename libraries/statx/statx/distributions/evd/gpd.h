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

#ifndef STATX_DISTRIBUTIONS_EVD_GPD_H_
#define STATX_DISTRIBUTIONS_EVD_GPD_H_

#include <cmath>
#include <vector>
#include <limits>
#include "statx/distributions/evd/common.h"

namespace statx {
namespace distributions {
namespace evd {
// Generalized Pareto Distribution (GPD)
// Evaluates the GPD density function for Maxima
// xi: shape parameter
// sigma: scale parameter
inline double gpdpdf(const double x, const double xi, const double sigma) {
  // Check for sigma > 0
  if (sigma <= 0) return std::numeric_limits<double>::infinity();
  if (x < 0) return 0.0;  // Defined x >= 0
  const double sigma_inv = 1.0/sigma;
  if (xi == 0.0) {  // Exponential case
    return sigma_inv*exp(-x*sigma_inv);
  }
  // Not exponential
  double arg = 1.0 + xi*x*sigma_inv;
  if (arg <= 0) return std::numeric_limits<double>::infinity();
  return sigma_inv*pow(arg, -1.0/xi - 1.0);
}

// Evaluates the GPD density function for minima
// xi: shape parameter
// sigma: scale parameter
inline double gpdcdf(const double x, const double xi, const double sigma) {
  // Check for sigma > 0
  if (sigma <= 0) return std::numeric_limits<double>::infinity();
  if (x < 0) return 0.0;  // Defined x >= 0
  const double sigma_inv = 1.0/sigma;
  if (xi == 0.0) {  // Exponential case
    return 1.0 - exp(-x*sigma_inv);
  }
  // Not exponential
  double arg = 1.0 + xi*x*sigma_inv;
  if (arg <= 0) return std::numeric_limits<double>::infinity();
  return 1.0 - pow(arg, -1.0/xi);
}

// Evaluates the pth-quantile
// xi: shape parameter
// sigma: scale parameter
inline double gpd_quantile(const double p,
                           const double xi,
                           const double sigma) {
  if (p < 0.0 || p > 1.0) return std::numeric_limits<double>::infinity();
  if (xi == 0) {  // Exponential
    return -sigma*log(1.0 - p);
  }
  // Non-exponential
  return sigma*(pow(1.0 - p, -xi) - 1.0)/xi;
}

/////////////////////////////
// Parameter estimation for the generalized Pareto distribution (GPD)
// xi: shape parameter
// sigma: scale parameter
bool gpdfit(const std::vector<double>& data,
              double* xi,
              double* sigma,
              FitType fit_type = MLE);
}  // evd
}  // distributions
}  // statx
#endif  // STATX_DISTRIBUTIONS_EVD_GPD_H_
