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

#ifndef STATX_DISTRIBUTIONS_EVD_GEV_H_
#define STATX_DISTRIBUTIONS_EVD_GEV_H_

#include <cmath>
#include <vector>
#include <limits>
#include "statx/distributions/evd/common.h"

namespace statx {
namespace distributions {
namespace evd {
// Generalized Extreme Value Distribution
// Evaluates the GEV density function for Maxima
// mu: location
// sigma: scale
// xi: shape
// pp212 castillo
inline double gevpdf(const double x,
                     const double mu,
                     const double sigma,
                     const double xi) {
  if (sigma <= 0) return std::numeric_limits<double>::infinity();
  double sigma_inv = 1.0/sigma;
  double arg = (x - mu)*sigma_inv;
  if (xi == 0.0) {  // Gumbel case
    double arg1 = exp(-arg);
    return sigma_inv*exp(-arg - arg1);
  }
  // Non-Gumbel case (Weibull/Frechet)
  double xi_inv = 1.0/xi;
  double temp1 = 1.0 + xi*arg;
  if (temp1 <= 0.0) return 0.0;
  double pow_factor = pow(temp1, -xi_inv - 1);
  double pow_factor2 = pow(temp1, -xi_inv);
  double factor1 = exp(-pow_factor2);
  return factor1*pow_factor*sigma_inv;
}

// Evalutes the GEV distribution function for Maxima
// mu: location
// sigma: scale
// xi: shape
inline double gevcdf(const double x,
                     const double mu,
                     const double sigma,
                     const double xi) {
  if (sigma <= 0) return std::numeric_limits<double>::infinity();
  double arg = (x - mu)/sigma;
  if (xi == 0.0) {  // Gumbel case
    return exp(-exp(-arg));
  }
  // Non-Gumbel case (Weibull/Frechet)
  double arg2 = 1 + xi*arg;
  if (arg2 <= 0.0) return 0.0;
  return exp(-pow(arg2, -1.0/xi));
}

// Evaluates the pth-quantile
// mu: location
// sigma: scale
// xi: shape
inline double gev_quantile(const double p,
                           const double mu,
                           const double sigma,
                           const double xi) {
  // p \in [0, 1]
  if (p < 0.0 || p > 1.0 || sigma <= 0) {
    return std::numeric_limits<double>::infinity();
  }
  double log_term = -log(p);
  if (xi == 0.0) {  // Gumbel case
    return mu - sigma*log(log_term);
  }
  return mu - sigma*(1.0 - pow(log_term, -xi))/xi;
}

// Find the parameters for the generalized extreme value (GEV) distribution by
// solving the MLE problem given the data. The distribution has three
// parameters:
//   mu: location
//   sigma: scale
//   xi: shape (determines the tail shape)
// When xi = 0, we have the case where the GEV is a Gumbel distribution.
bool gevfit(const std::vector<double>& data,
            double* mu,
            double* sigma,
            double* xi,
            FitType fit_type = MLE);
}  // evd
}  // distributions
}  // statx
#endif  // STATX_DISTRIBUTIONS_EVD_GEV_H_
