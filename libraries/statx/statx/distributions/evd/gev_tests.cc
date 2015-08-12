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

#include <glog/logging.h>
#include "gtest/gtest.h"
#include "statx/distributions/evd/gev.h"
#include "statx/distributions/evd/gev_mle.h"
#ifdef STATX_WITH_CERES
#include "statx/distributions/evd/gev_ceres.h"
#endif
#include "statx/utils/common_funcs.h"
#include "statx/utils/ecdf.h"

namespace statx {
namespace distributions {
namespace evd {
using std::vector;
using Eigen::Matrix;

// Matlab params: mu=10, sigma=2.5, xi=0.5
// R params ML estimates:
//        xi      sigma         mu
// 0.4841145  2.9570453 10.2058564
const vector<double> gev_data {
  1.182570e+01, 9.846248e+00, 1.344376e+01, 8.698939e+00, 6.962212e+00,
      1.158664e+01, 1.143223e+01, 1.147479e+01, 1.328162e+01, 2.322478e+01,
      2.586844e+01, 2.067551e+01, 1.204292e+01, 2.305931e+01, 8.413116e+00,
      9.819168e+00, 7.750681e+00, 1.135406e+01, 1.405570e+01, 8.825800e+00,
      1.041455e+01, 8.281900e+00, 9.577167e+00, 1.503409e+01, 1.134065e+01,
      1.475767e+01, 1.247069e+01, 1.987393e+01, 7.987021e+00, 8.792049e+00,
      1.034134e+01, 1.410769e+01, 1.985969e+01, 7.615216e+00, 8.550054e+00,
      1.039850e+01, 1.465032e+01, 1.122366e+01, 1.442014e+01, 8.754898e+00,
      1.294135e+01, 1.221529e+01, 7.239637e+00, 1.411656e+01, 5.988589e+01,
      8.488487e+00, 9.950656e+00, 8.886495e+00, 1.362177e+01, 8.805039e+00,
      4.845610e+01, 7.477011e+00, 7.784672e+00, 1.571509e+01, 1.786438e+01,
      1.163014e+01, 1.441198e+01, 1.117868e+01, 1.183843e+01, 8.884404e+00,
      1.103555e+01, 7.897587e+00, 7.945923e+00, 9.782476e+00, 1.236729e+01,
      1.978837e+01, 1.295134e+01, 1.313307e+01, 1.330143e+01, 1.557969e+01,
      1.276526e+01, 1.118251e+01, 9.759931e+00, 2.416089e+01, 9.225037e+00,
      1.110400e+01, 1.411845e+01, 9.272166e+00, 1.720868e+01, 1.135201e+01,
      2.129425e+01, 9.884607e+00, 2.433670e+01, 2.325427e+01, 1.185808e+01,
      8.307621e+00, 9.993185e+00, 9.405818e+00, 9.345539e+00, 9.103271e+00,
      6.941908e+00, 1.986355e+01, 1.062272e+01, 1.175766e+01, 9.655925e+00,
      6.954383e+01, 4.423742e+01, 3.131732e+01, 1.284641e+01, 1.383901e+01 };

// Matlab params: mu=5, sigma=2.5, xi=-0.4
const vector<double> gev_data2 {
  7.934858e+00, 8.772357e+00, 2.899042e+00, 8.858082e+00, 6.675548e+00,
      2.487359e+00, 4.354924e+00, 6.143097e+00, 9.467762e+00, 9.601250e+00,
      3.260413e+00, 9.715921e+00, 9.461948e+00, 5.760978e+00, 7.822001e+00,
      3.081604e+00, 5.356845e+00, 8.885563e+00, 7.760431e+00, 9.502272e+00,
      6.824108e+00, 1.134753e+00, 8.220760e+00, 9.113943e+00, 6.972406e+00,
      7.507758e+00, 7.404853e+00, 5.163404e+00, 6.822427e+00, 3.405271e+00,
      7.152216e+00, 9.965631e-01, 4.342702e+00, 1.454190e+00, 2.481040e+00,
      8.004974e+00, 7.077815e+00, 4.644002e+00, 9.348436e+00, 1.091083e+00,
      5.466153e+00, 5.092292e+00, 7.563474e+00, 7.783136e+00, 3.563525e+00,
      5.788422e+00, 5.509854e+00, 6.763969e+00, 7.174389e+00, 7.486065e+00,
      4.335719e+00, 6.978702e+00, 6.819997e+00, 3.314695e+00, 2.794840e+00,
      5.842088e+00, 9.506712e+00, 4.810157e+00, 6.380952e+00, 3.905490e+00,
      7.461888e+00, 4.169313e+00, 5.889365e+00, 7.105896e+00, 8.614020e+00,
      9.498735e+00, 6.145165e+00, 3.042826e+00, 3.167431e+00, 4.188873e+00,
      8.148301e+00, 4.162700e+00, 7.931374e+00, 4.074043e+00, 9.051784e+00,
      4.877146e+00, 3.657353e+00, 4.136567e+00, 6.572923e+00, 5.685167e+00,
      4.888760e+00, 8.065363e+00, 6.380929e+00, 6.160687e+00, 8.902751e+00,
      4.411402e+00, 7.503915e+00, 7.479282e+00, 5.084833e+00, 6.272737e+00,
      2.120298e+00, 1.655662e+00, 6.043519e+00, 7.663021e+00, 9.114176e+00,
      2.935956e+00, 6.278945e+00, 5.660638e+00, -8.658902e-02, 4.787194e+00 };

const vector<double> gevpdf_data {
  0, 0, 2.300535e-268,  4.268253e-118, 1.082125e-65, 1.488030e-41, 1.603900e-28,
      1.013431e-20, 1.059429e-15, 2.704530e-12, 6.943972e-10, 3.997288e-08,
      8.348174e-07, 8.563777e-06, 5.261571e-05, 2.214124e-04, 7.005724e-04,
      1.781355e-03, 3.820453e-03, 7.163573e-03, 1.206534e-02, 1.863482e-02,
      2.681917e-02, 3.642127e-02, 4.713870e-02, 5.861004e-02, 7.045816e-02,
      8.232441e-02, 9.389194e-02, 1.048986e-01, 1.151417e-01, 1.244763e-01,
      1.328100e-01, 1.400958e-01, 1.463238e-01, 1.515132e-01, 1.557052e-01,
      1.589563e-01, 1.613334e-01, 1.629092e-01, 1.637589e-01, 1.639575e-01,
      1.635782e-01, 1.626907e-01, 1.613607e-01, 1.596491e-01, 1.576118e-01,
      1.552997e-01, 1.527589e-01, 1.500307e-01, 1.471518e-01, 1.441548e-01,
      1.410686e-01, 1.379183e-01, 1.347259e-01, 1.315106e-01, 1.282889e-01,
      1.250748e-01, 1.218806e-01, 1.187163e-01, 1.155907e-01, 1.125108e-01,
      1.094827e-01, 1.065110e-01, 1.035995e-01, 1.007514e-01, 9.796873e-02,
      9.525322e-02, 9.260591e-02, 9.002742e-02, 8.751794e-02, 8.507733e-02,
      8.270517e-02, 8.040080e-02, 7.816333e-02, 7.599175e-02, 7.388487e-02,
      7.184141e-02, 6.986000e-02, 6.793920e-02, 6.607752e-02, 6.427345e-02,
      6.252543e-02, 6.083190e-02, 5.919131e-02, 5.760210e-02, 5.606272e-02,
      5.457165e-02, 5.312739e-02, 5.172844e-02, 5.037337e-02, 4.906074e-02,
      4.778918e-02, 4.655732e-02, 4.536386e-02, 4.420751e-02, 4.308701e-02,
      4.200117e-02, 4.094882e-02, 3.992881e-02, 3.894004e-02 };

const vector<double> gevcdf_data {
  0, 0, 3.680856e-272, 2.304857e-121, 1.385119e-68, 3.720076e-44, 6.928847e-31,
      6.952136e-23, 1.084855e-17, 3.943205e-14, 1.388794e-11, 1.064078e-09,
      2.885129e-08, 3.762924e-07, 2.887550e-06, 1.494534e-05, 5.739089e-05,
      1.750360e-04, 4.456176e-04, 9.826989e-04, 1.930454e-03, 3.451542e-03,
      5.711411e-03, 8.862752e-03, 1.303291e-02, 1.831564e-02, 2.476745e-02,
      3.240783e-02, 4.122232e-02, 5.116745e-02, 6.217652e-02, 7.416545e-02,
      8.703837e-02, 1.006925e-01, 1.150222e-01, 1.299226e-01, 1.452916e-01,
      1.610323e-01, 1.770537e-01, 1.932722e-01, 2.096114e-01, 2.260023e-01,
      2.423836e-01, 2.587010e-01, 2.749070e-01, 2.909605e-01, 3.068260e-01,
      3.224737e-01, 3.378783e-01, 3.530192e-01, 3.678794e-01, 3.824456e-01,
      3.967074e-01, 4.106572e-01, 4.242897e-01, 4.376016e-01, 4.505916e-01,
      4.632597e-01, 4.756072e-01, 4.876368e-01, 4.993518e-01, 5.107565e-01,
      5.218557e-01, 5.326549e-01, 5.431599e-01, 5.533769e-01, 5.633123e-01,
      5.729729e-01, 5.823653e-01, 5.914963e-01, 6.003730e-01, 6.090022e-01,
      6.173908e-01, 6.255455e-01, 6.334732e-01, 6.411804e-01, 6.486737e-01,
      6.559595e-01, 6.630440e-01, 6.699335e-01, 6.766338e-01, 6.831509e-01,
      6.894904e-01, 6.956578e-01, 7.016585e-01, 7.074978e-01, 7.131806e-01,
      7.187120e-01, 7.240965e-01, 7.293389e-01, 7.344437e-01, 7.394150e-01,
      7.442572e-01, 7.489742e-01, 7.535699e-01, 7.580482e-01, 7.624126e-01,
      7.666668e-01, 7.708140e-01, 7.748576e-01, 7.788008e-01 };

// mu=5, sigma=2.5, xi=0.5
const vector<double> gevinv_data {
  0, 3.295051e+00, 3.941240e+00, 4.556818e+00, 5.223401e+00, 6.005612e+00,
      6.995740e+00, 8.372086e+00, 1.058468e+01, 1.540391e+01 };

TEST(GEV, MLE1) {
  const int n = gev_data.size() + 4;
  GEVMLEObjective mle(gev_data);
  Matrix<double, Dynamic, 1> x(n);
  x.setZero();

  // True Values
  x(0) = 10;  // mu (Location parameter)
  x(1) = 2.5;  // sigma (Scale parameter)
  x(2) = 0.5;  // xi (Shape parameter)
  // The lagrange multipliers are set to zero.
  const double mle_matlab = 296.1117;
  double mle_val = mle(x);
  VLOG(1) << "MLE=" << mle_val << " MLE (Matlab): " << mle_matlab;
  ASSERT_NEAR(mle_matlab, mle_val, 1e-3);
}

TEST(GEV, MLE2) {
  const int n = gev_data2.size() + 4;
  GEVMLEObjective mle(gev_data2);
  Matrix<double, Dynamic, 1> x(n);
  x.setZero();

  // True Values
  x(0) = 5.0;  // mu (Location parameter)
  x(1) = 2.5;  // sigma (Scale parameter)
  x(2) = -0.4;  // xi (Shape parameter)
  // The lagrange multipliers are set to zero.
  const double mle_matlab = 224.0336;
  double mle_val = mle(x);
  VLOG(1) << "MLE=" << mle_val << " MLE (Matlab): " << mle_matlab;
  ASSERT_NEAR(mle_matlab, mle_val, 1e-3);
}

TEST(GEV, MLE_Penalties) {
  const int n = 3;
  GEVMLEObjective mle(gev_data);
  Matrix<double, Dynamic, 1> x(n);
  // True Values
  x(0) = 10;  // mu (Location parameter)
  x(1) = 1e-4;
  x(2) = 0.5;  // xi (Shape parameter)
  double mle_val = mle(x);
  VLOG(1) << "MLE val: " << mle_val;
  ASSERT_GT(mle_val, 1e+10);

  // Bad Values that should add some cost
  // Violating 1 + xi*(z - mu)/sigma
  x(0) = 10;  // mu (Location parameter)
  x(1) = 0.01;  // sigma (Scale parameter)
  x(2) = -0.5;  // xi (Shape parameter)
  mle_val = mle(x);
  VLOG(1) << "MLE val: " << mle_val;
  ASSERT_GT(mle_val, 1e+10);

  // Bad Values that should add some cost
  // Violating sigma <= 0
  x(0) = 10;  // mu (Location parameter)
  x(1) = -1.0;  // sigma (Scale parameter)
  x(2) = 0.5;  // xi (Shape parameter)
  mle_val = mle(x);
  VLOG(1) << "MLE val: " << mle_val;
  ASSERT_GT(mle_val, 0.0);
}

TEST(GEV, GradientFunctor) {
  // int n = gev_data.size() + 4;
  const int n = 3;
  GEVMLEGradientFunctor gradient(gev_data);
  Matrix<double, Dynamic, 1> x(n), g(n);

  // ML estimates (MATLAB): xi=0.4841 sigma=2.9570  mu=10.2059
  // Should be really close to a *zero* norm.
  x(0) = 10.2059;  // mu (Location parameter)
  x(1) = 2.9570;  // sigma (Scale parameter)
  x(2) = 0.4841;  // xi (Shape parameter)
  gradient(x, &g);
  VLOG(1) << "Gradient: " << g.transpose();
  double norm = g.norm();
  ASSERT_LT(norm, sqrt(n));

  // ML estimates (MATLAB): xi=-0.5563 sigma=2.5907 mu=5.4223
  GEVMLEGradientFunctor gradient2(gev_data2);
  x.setConstant(0.0);
  g.setConstant(0.0);
  // True Values
  x(0) = 5.4223;  // mu (Location parameter)
  x(1) = 2.5907;  // sigma (Scale parameter)
  x(2) = -0.5563;  // xi (Shape parameter)
  gradient2(x, &g);
  VLOG(1) << "Gradient: " << g.transpose();
  norm = g.norm();
  ASSERT_LT(norm, sqrt(n));
}

TEST(GEV, PDF) {
  // Params
  const double mu = 5.0;
  const double sigma = 2.5;
  const double xi = 0.5;
  // Domain
  const int nsamples = 100;
  const double x_final = 10.0;
  const double x_start = 0.0;
  const double dx = (x_final - x_start)/nsamples;
  // Testing
  for (int i = 0; i < nsamples; i++) {
    double x = i*dx;
    double y_gt = gevpdf_data[i];
    double y = gevpdf(x, mu, sigma, xi);
    ASSERT_NEAR(y_gt, y, 1e-2);
  }
  // Test for last samples (x_final)
  double y_gt = gevpdf_data[nsamples];
  double y = gevpdf(x_final, mu, sigma, xi);
  ASSERT_NEAR(y_gt, y, 1e-2);
}

TEST(GEV, CDF) {
  // Params
  const double mu = 5.0;
  const double sigma = 2.5;
  const double xi = 0.5;
  // Domain
  const int nsamples = 100;
  const double x_final = 10.0;
  const double x_start = 0.0;
  const double dx = (x_final - x_start)/nsamples;
  // Testing
  for (int i = 0; i < nsamples; i++) {
    double x = i*dx;
    double y_gt = gevcdf_data[i];
    double y = gevcdf(x, mu, sigma, xi);
    ASSERT_NEAR(y_gt, y, 1e-2);
  }
  // Test for last samples (x_final)
  double y_gt = gevcdf_data[nsamples];
  double y = gevcdf(x_final, mu, sigma, xi);
  ASSERT_NEAR(y_gt, y, 1e-2);
}

TEST(GEV, FitMLE1) {
  const double mu_gt = 10.0;
  const double sigma_gt = 2.5;
  const double xi_gt = 0.5;
  double mu = 10.0;
  double sigma = 2.5;
  double xi = 0.5;
  ASSERT_TRUE(gevfit(gev_data, &mu, &sigma, &xi));
  VLOG(1) << "mu=" << mu << " sigma=" << sigma << " xi=" << xi;
  ASSERT_NEAR(mu_gt, mu, 1.0);
  ASSERT_NEAR(sigma_gt, sigma, 1.0);
  ASSERT_NEAR(xi_gt, xi, 1.0);
}

TEST(GEV, FitMLE2) {
  const double mu_gt = 5.0;
  const double sigma_gt = 2.5;
  const double xi_gt = -0.4;
  double mu = 10.0;
  double sigma = 2.5;
  double xi = 0.5;
  ASSERT_TRUE(gevfit(gev_data2, &mu, &sigma, &xi));
  VLOG(1) << "mu=" << mu << " sigma=" << sigma << " xi=" << xi;
  ASSERT_NEAR(mu_gt, mu, 1.0);
  ASSERT_NEAR(sigma_gt, sigma, 1.0);
  ASSERT_NEAR(xi_gt, xi, 1.0);
}

TEST(GEV, Quantile) {
  // Params
  const double mu = 5.0;
  const double sigma = 2.5;
  const double xi = 0.5;

  double p = 0.0;
  for (int i = 0; i < gevinv_data.size(); i++) {
    const double z = gevinv_data[i];
    double q = gev_quantile(p, mu, sigma, xi);
    ASSERT_NEAR(q, z, 0.1);
    p += 0.1;
  }
}

#ifdef STATX_WITH_CERES
//////////////////////////////////////////
// Testing GEV fit w/ CERES
TEST(GEV, CERES_GEVCostFunction) {
  // Calculate the ECDF from the data
  std::vector<double> fx, x;
  statx::utils::ecdf(gev_data, &fx, &x);
  const double var = statx::utils::stddev(gev_data);
  const double sigma0 = sqrt(6*var*var)/M_PI;  // sigma
  const double mu0 = statx::utils::mean(gev_data) -0.57722*sigma0;  // mu
  const double xi0 = 0.1;  // according to EVIR's package

  double mu = mu0;
  double sigma = sigma0;
  double xi = xi0;
  double residual;
  double* params = new double[3];
  params[0] = mu;
  params[1] = sigma;
  params[2] = xi;
  double* parameters[1];
  parameters[0] = params;
  double* jacobians[1];
  double* jacobian = new double[3];
  jacobians[0] = jacobian;
  double residuals;

  for (int i = 0; i < fx.size() - 1; i++) {
    GEVCostFunctionAnalytic cost(x[i], fx[i]);
    ASSERT_TRUE(cost.Evaluate(&parameters[0], &residuals, &jacobians[0]));
    ASSERT_TRUE(std::isfinite(residual));
  }

  delete [] params;
  params = nullptr;
  delete [] jacobian;
  jacobian = nullptr;
}

TEST(GEV, FitQuantileLeastSquares_Case1) {
  double mu = 10.0;
  double sigma = 2.5;
  double xi = 0.5;
  ASSERT_TRUE(gevfit(gev_data, &mu, &sigma, &xi, QUANTILE_NLS));
  VLOG(1) << "mu=" << mu << " sigma=" << sigma << " xi=" << xi;
  vector<double> fx, x;
  statx::utils::ecdf(gev_data, &fx, &x);
  double mse = 0.0;

  for (int i = 0; i < x.size(); i++) {
    double y = gevcdf(x[i], mu, sigma, xi);
    double error = y - fx[i];
    mse += error*error;
  }
  mse = sqrt(mse)/x.size();
  VLOG(1) << "MSE: " << mse;
  ASSERT_LT(mse, 1.0);

  const double mu_gt = 10.0;
  const double sigma_gt = 2.5;
  const double xi_gt = 0.5;
  ASSERT_NEAR(mu_gt, mu, 1.0);
  ASSERT_NEAR(sigma_gt, sigma, 1.0);
  ASSERT_NEAR(xi_gt, xi, 1.0);
}

TEST(GEV, FitQuantileLeastSquares_Case2) {
  double mu = 5;
  double sigma = 2.5;
  double xi = -0.3;

  ASSERT_TRUE(gevfit(gev_data2, &mu, &sigma, &xi, QUANTILE_NLS));
  VLOG(1) << "mu=" << mu << " sigma=" << sigma << " xi=" << xi;

  vector<double> fx, x;
  statx::utils::ecdf(gev_data2, &fx, &x);
  double mse = 0.0;

  for (int i = 0; i < x.size(); i++) {
    double y = gevcdf(x[i], mu, sigma, xi);
    double error = y - fx[i];
    mse += error*error;
  }
  mse = sqrt(mse)/x.size();
  VLOG(1) << "MSE: " << mse;
  ASSERT_LT(mse, 1.0);

  const double mu_gt = 5.0;
  const double sigma_gt = 2.5;
  const double xi_gt = -0.4;
  ASSERT_NEAR(mu_gt, mu, 1.0);
  ASSERT_NEAR(sigma_gt, sigma, 1.0);
  ASSERT_NEAR(xi_gt, xi, 1.0);
}
#endif
}  // evd
}  // distributions
}  // statx
