// Copyright (C) 2013  Victor Fragoso <vfragoso@cs.ucsb.edu>
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

#include "statx/utils/ecdf.h"
#include <cmath>
#include <random>
#include <glog/logging.h>
#include "gtest/gtest.h"

namespace statx {
namespace utils {
TEST(ECDF, NO_DATA) {
  // Calculating the ECDF
  std::vector<double> fx, x, y;
  ecdf(y, &fx, &x);
  ASSERT_EQ(fx.size(), 0);

  ecdf(y, &fx, &x);
  ASSERT_EQ(x.size(), 0);
}

TEST(ECDF, RAYLEIGH_CDF) {
  // Generate data for Rayleigh
  std::random_device rd;
  std::mt19937 rng(rd());
  double lambda = 1.0;  // Exponential parameter
  double sigma = 1.0;  // Rayleigh parameter
  std::exponential_distribution<> exp_dist(lambda);
  const int nsamples = 1000;

  std::vector<double> data(nsamples);
  for (int i = 0; i < nsamples; i++) {
    // Rayleigh = sqrt(2*exponential*sigma^2*lambda)
    data[i] = sqrt(2*exp_dist(rng)*sigma*sigma*lambda);
  }

  // Build ECDF
  std::vector<double> fx, x;
  ecdf(data, &fx, &x);

  for (unsigned int i = 0; i < fx.size(); i++) {
    double F = 1.0 - exp(-x[i]*x[i]/2);  // CDF closed form
    double err = F - fx[i];
    ASSERT_GT(0.9, fabs(err));
  }
}

TEST(ECDF, EXPONENTIAL_CDF) {
  // Generate data for Exponential
  std::random_device rd;
  std::mt19937 rng(rd());
  double lambda = 1.0;  // Exponential parameter
  std::exponential_distribution<> exp_dist(lambda);
  const int nsamples = 1000;

  std::vector<double> data(nsamples);
  for (int i = 0; i < nsamples; i++) {
    data[i] = exp_dist(rng);
  }

  // Build ECDF
  std::vector<double> fx, x;
  ecdf(data, &fx, &x);

  for (unsigned int i = 0; i < fx.size(); i++) {
    double F = 1.0 - exp(-lambda*x[i]);  // CDF closed form
    double err = F - fx[i];
    ASSERT_GT(0.9, fabs(err));
  }
}
}  // utils
}  // statx
