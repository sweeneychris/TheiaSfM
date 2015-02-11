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

#include "statx/distributions/rayleigh.h"
#include <glog/logging.h>
#include "gtest/gtest.h"

namespace statx {
namespace distributions {
using std::vector;
// Rayleigh data generated w/ z = raylrnd(5, 1, 100);
const vector<double> rayl_data {
  4.988279e+00, 1.018790e+01, 1.130532e+01, 5.098663e+00, 2.200773e+00,
      7.194513e+00, 3.271346e+00, 4.074461e+00, 1.983399e+01, 1.388116e+01,
      1.264393e+01, 1.574458e+01, 7.682996e+00, 5.370039e+00, 5.988092e+00,
      1.197909e+00, 7.210251e+00, 1.231293e+01, 7.114181e+00, 9.310377e+00,
      1.492265e+01, 7.312679e+00, 7.771736e+00, 9.717783e+00, 3.386066e+00,
      5.349835e+00, 6.585778e+00, 2.057238e+00, 3.803028e+00, 1.098836e+01,
      4.781361e+00, 7.060568e+00, 9.525745e+00, 4.778406e+00, 1.478878e+01,
      7.193849e+00, 6.863756e+00, 6.784459e+00, 7.071675e+00, 8.687245e+00,
      5.238669e-01, 1.781420e+00, 8.895415e+00, 2.118236e+00, 5.998307e+00,
      4.898337e+00, 5.840483e+00, 4.119161e+00, 1.140809e+01, 7.348853e+00,
      5.040731e+00, 3.995480e-01, 6.073074e+00, 6.850089e+00, 5.093542e+00,
      7.692045e+00, 5.251045e+00, 7.007435e+00, 1.592250e+00, 6.315536e+00,
      5.639930e+00, 4.242755e+00, 6.244945e+00, 1.377681e+01, 1.131952e+01,
      1.596575e+00, 9.753446e+00, 5.700980e+00, 5.380799e+00, 1.240059e+01,
      7.336896e+00, 1.223509e+01, 7.309083e+00, 4.745677e+00, 4.295789e+00,
      7.366092e+00, 7.141651e+00, 2.608040e+00, 2.540589e+00, 4.420176e+00,
      9.834446e+00, 2.702191e+00, 4.357284e+00, 8.113165e+00, 4.880407e+00,
      6.254988e+00, 6.975831e+00, 6.420003e+00, 3.472400e+00, 5.839088e+00,
      5.822171e+00, 1.207800e+00, 3.916392e+00, 1.292797e+01, 3.344301e+00,
      4.235188e+00, 7.646017e+00, 9.943071e+00, 2.432442e+00, 9.502283e+00 };

// PDF calculated w/ d = raylpdf(x, sigma); for
// x = 0:0.1:10 and sigma = 5
const vector<double> rayl_pdf_data {
  0, 3.999200e-03, 7.993603e-03, 1.197842e-02, 1.594888e-02, 1.990025e-02,
      2.382782e-02, 2.772694e-02, 3.159301e-02, 3.542150e-02, 3.920795e-02,
      4.294798e-02, 4.663732e-02, 5.027177e-02, 5.384727e-02, 5.735985e-02,
      6.080567e-02, 6.418103e-02, 6.748235e-02, 7.070621e-02, 7.384931e-02,
      7.690853e-02, 7.988090e-02, 8.276362e-02, 8.555404e-02, 8.824969e-02,
      9.084828e-02, 9.334770e-02, 9.574600e-02, 9.804143e-02, 1.002324e-01,
      1.023176e-01, 1.042957e-01, 1.061658e-01, 1.079270e-01, 1.095786e-01,
      1.111203e-01, 1.125517e-01, 1.138726e-01, 1.150833e-01, 1.161838e-01,
      1.171748e-01, 1.180566e-01, 1.188301e-01, 1.194961e-01, 1.200558e-01,
      1.205104e-01, 1.208611e-01, 1.211095e-01, 1.212573e-01, 1.213061e-01,
      1.212579e-01, 1.211147e-01, 1.208785e-01, 1.205517e-01, 1.201364e-01,
      1.196351e-01, 1.190502e-01, 1.183844e-01, 1.176403e-01, 1.168205e-01,
      1.159279e-01, 1.149651e-01, 1.139351e-01, 1.128407e-01, 1.116849e-01,
      1.104706e-01, 1.092006e-01, 1.078781e-01, 1.065060e-01, 1.050871e-01,
      1.036245e-01, 1.021212e-01, 1.005800e-01, 9.900392e-02, 9.739574e-02,
      9.575834e-02, 9.409452e-02, 9.240704e-02, 9.069863e-02, 8.897194e-02,
      8.722958e-02, 8.547412e-02, 8.370803e-02, 8.193376e-02, 8.015367e-02,
      7.837004e-02, 7.658510e-02, 7.480099e-02, 7.301981e-02, 7.124353e-02,
      6.947409e-02, 6.771332e-02, 6.596299e-02, 6.422478e-02, 6.250029e-02,
      6.079105e-02, 5.909849e-02, 5.742397e-02, 5.576878e-02, 5.413411e-02 };

// PDF calculated w/ d = raylcdf(x, sigma); for
// x = 0:0.1:10 and sigma = 5
const vector<double> rayl_cdf_data {
  0, 1.999800e-04, 7.996801e-04, 1.798381e-03, 3.194885e-03, 4.987521e-03,
      7.174142e-03, 9.752136e-03, 1.271843e-02, 1.606949e-02, 1.980133e-02,
      2.390953e-02, 2.838923e-02, 3.323516e-02, 3.844162e-02, 4.400252e-02,
      4.991137e-02, 5.616130e-02, 6.274510e-02, 6.965519e-02, 7.688365e-02,
      8.442226e-02, 9.226246e-02, 1.003954e-01, 1.088121e-01, 1.175031e-01,
      1.264588e-01, 1.356694e-01, 1.451250e-01, 1.548152e-01, 1.647298e-01,
      1.748582e-01, 1.851897e-01, 1.957137e-01, 2.064193e-01, 2.172955e-01,
      2.283313e-01, 2.395158e-01, 2.508380e-01, 2.622867e-01, 2.738510e-01,
      2.855198e-01, 2.972823e-01, 3.091275e-01, 3.210447e-01, 3.330232e-01,
      3.450524e-01, 3.571218e-01, 3.692212e-01, 3.813404e-01, 3.934693e-01,
      4.055983e-01, 4.177178e-01, 4.298182e-01, 4.418904e-01, 4.539256e-01,
      4.659149e-01, 4.778498e-01, 4.897222e-01, 5.015241e-01, 5.132477e-01,
      5.248858e-01, 5.364310e-01, 5.478765e-01, 5.592159e-01, 5.704426e-01,
      5.815509e-01, 5.925349e-01, 6.033893e-01, 6.141089e-01, 6.246889e-01,
      6.351248e-01, 6.454125e-01, 6.555478e-01, 6.655273e-01, 6.753475e-01,
      6.850055e-01, 6.944983e-01, 7.038236e-01, 7.129790e-01, 7.219627e-01,
      7.307729e-01, 7.394082e-01, 7.478674e-01, 7.561495e-01, 7.642539e-01,
      7.721801e-01, 7.799279e-01, 7.874972e-01, 7.948882e-01, 8.021013e-01,
      8.091371e-01, 8.159964e-01, 8.226801e-01, 8.291894e-01, 8.355255e-01,
      8.416900e-01, 8.476843e-01, 8.535103e-01, 8.591697e-01, 8.646647e-01 };

TEST(Rayleigh, PDF) {
  double x = 0.0;
  const double sigma = 5.0;
  for (int i = 0; i < rayl_pdf_data.size(); i++) {
    const double d = raylpdf(x, sigma);
    ASSERT_NEAR(rayl_pdf_data[i], d, 1e-3);
    x += 0.1;
  }
}

TEST(Rayleigh, CDF) {
  double x = 0.0;
  const double sigma = 5.0;
  for (int i = 0; i < rayl_cdf_data.size(); i++) {
    const double d = raylcdf(x, sigma);
    ASSERT_NEAR(rayl_cdf_data[i], d, 1e-3);
    x += 0.1;
  }
}

TEST(Rayleigh, FitMLE) {
  const double sigma_gt = 5.0;
  const double sigma = raylfit(rayl_data);
  ASSERT_NEAR(sigma_gt, sigma, 1.0);
}
}  // distributions
}  // statx
