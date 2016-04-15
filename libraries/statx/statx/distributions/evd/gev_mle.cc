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

#include "statx/distributions/evd/gev_mle.h"
#include <cmath>
#include <vector>
#include <limits>
#include <glog/logging.h>
#include "statx/utils/common_funcs.h"
#include <optimo/solvers/bfgs.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif


namespace statx {
namespace distributions {
namespace evd {

using Eigen::Matrix;
using Eigen::Dynamic;
using std::vector;

// MLE objective function for the generalized extreme value distribution (GEV)
// The distribution has three parameters:
// -\inf < mu < \inf : location
// 0 < sigma : scale
// -\inf < xi < \inf: shape
// When xi = 0 => the GEV becomes the Gumbel family of distributions
// We implement the negative log-likelihood function.
// The problem we try to solve given m observations is:
//        minimize -l(z; mu, sigma, xi)
//
// The minimizer is the vector holding our ML estimates for the GEV distribution
// given the data.
// There are two log-likelihood functions that need to be optimized. The
// selection of which function to optimize depends on the xi. If xi = 0, then
// we solve for Gumbel family of distributions. Otherwise, use the general
// formulation of the log-likelihood.
// Note: It is well known that the MLE estimator for the GEV distribution can
// have a bad behavior. These are the identified cases:
// 1. When xi > -0.5, the  ML properties are regular.
// 2. When -1 < xi < -0.5, ML are often obtainable but lack the asymptotic
// ML properties.
// 3. When xi < -1, the ML are unlikely to be obtain

// Implementation for xi != 0
// Problem object
class GEVMLEProblem : public optimo::ProblemLS<double> {
 public:
  // Constructor
  GEVMLEProblem(const GEVMLEObjective& obj,
                const GEVMLEGradientFunctor& g,
                const GEVMLEHessianFunctor& h) :
      optimo::ProblemLS<double>(obj, g, h) { }

  // Destructor
  virtual ~GEVMLEProblem(void) { }
};

// Negative log-likelihood for GEV when xi != 0
// l(z; \mu, \sigma, \xi) = m*log(\sigma) +
//  (1 + \frac{1}{xi})\sum_{i=1}^m \log{1 + \xi \frac{z_i - \mu}{\sigma}} +
//  {1 + \xi \frac{z_i - \mu}{\sigma}}^{-\frac{1}{\xi}}
// The Lagrangian for this MLE problem is:
// L = l(z; \mu, \sigma, \xi) -\lambda_1 \log{\sigma} -
//   \sum_{i=2}^{m+1} \lambda_i \log{1.0 + \xi(z_{i-1} - \mu)/\sigma}
double GEVMLEObjective::operator()(const Matrix<double, Dynamic, 1>& x) const {
  const int m = data_.size();
  const double& mu = x(0);  // Location
  const double& sigma = x(1);  // Scale
  const double& xi = x(2);  // Shape
  const double lambda_inv = 1.0/lambda_;
  const double xi_term = 1.0 + 1.0/xi;
  double acc1 = 0.0;
  double acc2 = 0.0;
  double acc3 = 0.0;
  if (sigma <= 0.0) return std::numeric_limits<double>::max();
  for (double z : data_) {
    const double z_minus_mu = z - mu;
    const double z_minus_mu_sigma = z_minus_mu/sigma;
    const double log_arg = 1.0 + xi*z_minus_mu_sigma;
    // TODO(vfragoso): Should we keep the original condition?
    // Condition: 1 + xi*(z - mu)/sigma > 0
    // const double log_barrier_term = log(sigma + xi*z_minus_mu);
    if (log_arg < 0) return std::numeric_limits<double>::max();
    const double log_log_arg = log(log_arg);
    acc1 += log_log_arg;
    acc2 += pow(log_arg, -1.0/xi);
    // TODO(vfragoso): Document penalty
    acc3 += (log_arg <= lambda_inv) ? log(lambda_*log_arg) : 0.0;
  }

  // TODO(vfragoso): Document penalty
  const double penalty1 =
      (sigma <= lambda_inv) ? - alpha_*log(lambda_*sigma) : 0.0;
  return m*log(sigma) + xi_term*acc1 + acc2 + penalty1 - beta_*acc3;
}

// TODO(vfragoso): Rewrite to give more space and check that the equations
// are right.
// Gradient of the Negative log-likelihood from above is:
// g = [dl/d\mu dl/d\sigma dl/d\xi]'
// dl/d\mu = -(1 + \frac{1}{\xi})\sum_{i=1}^m\frac{1}{\sigma + \xi(z_i - \mu)} +
// \frac{1}{\sigma\xi}\sum_{i=1}^m (1 + \frac{1}{\xi})^{\frac{-1}{\xi} - 1}
// dl/d\sigma = -\frac{m}{\sigma} +
// \frac{1 + 1/\xi}{\sigma}\sum_{i=1}^m\frac{-z_i+\mu)}{\sigma + \xi(\z_i-\mu)}+
// \frac{1}{\xi\sigma^2}\sum_{i=1}^m(z_i-\mu)(1+\xi(\zi-\mu)/\sigma)^{-1/\xi-1}
// dl/d\xi = -\frac{1}{\xi^2}\sum_{i=1}^m\log{1+\xi(z_i-\mu)/\sigma} +
// (1+\frac{1}{\xi})\sum_{i=1}^m\frac{z_i-\mu}{\sigma + \xi(z_i-\mu)} +
// \sum_{i=1}^m f'(\xi)
// f'(\xi) = f(\xi)(1 + \xi(z_i-\mu)/\sigma)^{\frac{-1}{\xi}}
// f(\xi) = \frac{1}{\xi^2}\lig{1+\xi(\z_i-\mu)/\sigma} -
// \frac{z_i-\mu}{\xi\sigma(\sigma + \xi(z_i-\mu))}
void GEVMLEGradientFunctor::operator()(const Matrix<double, Dynamic, 1>& x,
                                       Matrix<double, Dynamic, 1>* g) const {
  const int m = data_.size();
  const double& mu = x(0);
  const double& sigma = x(1);
  const double& xi = x(2);
  const double sigma_inv = 1.0/sigma;
  const double xi_inv = 1.0/xi;
  const double sigma_sqrd_inv = sigma_inv*sigma_inv;
  const double xi_sqrd_inv = xi_inv*xi_inv;
  const double lambda_inv = 1.0/lambda_;
  const double sigma_sqrd = sigma*sigma;
  double mu_acc1 = 0.0;
  double mu_acc2 = 0.0;
  double sigma_acc1 = 0.0;
  double sigma_acc2 = 0.0;
  double xi_acc1 = 0.0;
  double xi_acc2 = 0.0;
  double xi_acc3 = 0.0;
  double acc1 = 0.0;
  double acc2 = 0.0;
  double acc3 = 0.0;

  for (double z : data_) {
    const double z_minus_mu = z - mu;
    const double xi_z_minus_mu = xi*z_minus_mu;
    const double xi_z_minus_mu_inv = 1.0 / xi_z_minus_mu;
    const double sigma_xi_z_minus_mu = sigma + xi_z_minus_mu;
    const double sigma_xi_z_minus_mu_inv = 1.0 / sigma_xi_z_minus_mu;
    // discriminant > 0 !
    const double discriminant = 1.0 + xi_z_minus_mu*sigma_inv;
    const double pow_term = pow(discriminant, -xi_inv - 1.0);
    const double pow_term2 = pow(discriminant, -xi_inv);
    const double log_term = log(discriminant);

    // Mu terms
    mu_acc1 += sigma_xi_z_minus_mu_inv;
    mu_acc2 += pow_term;

    // Sigma terms
    sigma_acc1 += z_minus_mu*sigma_xi_z_minus_mu_inv;
    sigma_acc2 += z_minus_mu*pow_term;

    // xi  terms
    xi_acc1 += log_term;
    xi_acc2 += z_minus_mu*sigma_xi_z_minus_mu_inv;
    xi_acc3 += pow_term2*(xi_inv*log_term -
                          z_minus_mu*sigma_xi_z_minus_mu_inv);

    // Terms for the constraints
    bool penalty = discriminant <= lambda_inv;
    const double lambda_d_inv = 1.0/lambda_*discriminant;
    acc1 += penalty ? lambda_d_inv : 0.0;
    acc2 += penalty ? z_minus_mu*lambda_d_inv : 0.0;
    acc3 += penalty ? z_minus_mu*lambda_d_inv : 0.0;
  }
  const double xi_plus_one = xi + 1.0;
  // Gradient from the neg-log-likelihood
  (*g)(0) = -xi_plus_one*mu_acc1 + sigma_inv*mu_acc2;
  (*g)(1) = m*sigma_inv - xi_plus_one*sigma_inv*sigma_acc1 +
      sigma_sqrd_inv*sigma_acc2;
  (*g)(2) = -xi_sqrd_inv*xi_acc1 + (1 + xi_inv)*xi_acc2 + xi_inv*xi_acc3;

  // Adding terms coming from the constraints (penalty terms)
  (*g)(0) += beta_*xi*acc1/sigma;
  (*g)(1) += beta_*acc2/sigma_sqrd;
  (*g)(2) -= beta_*acc3/sigma;
}

// Hessian (Using an approximated hessian from the gradients). Thus Newton
// becomes the BFGS method.
void
GEVMLEHessianFunctor::operator()(const Matrix<double, Dynamic, 1>& x,
                                 Matrix<double, Dynamic, Dynamic>* h) const {
}


bool gevfit_mle(const vector<double>& data,
                double* mu,
                double* sigma,
                double* xi) {
  bool exit_flag = false;
  // const int n = data.size() + 4;
  const int n = 3;
  Matrix<double, Dynamic, 1> x(n);
  static const double delta = 0.1;
  const double mean = statx::utils::mean(data);
  const double var = statx::utils::stddev(data, mean);
  const double sigma0 = sqrt(6*var*var)/M_PI;  // sigma
  const double mu0 = mean - 0.57722*sigma0;  // mu
  const double xi0 = 0.1;  // according to EVIR's package

  x.setConstant(0.0);
  x(0) = mu0;
  x(1) = sigma0;
  x(2) = xi0;

  // Solve using BFGS
  GEVMLEObjective mle(data);
  GEVMLEGradientFunctor gradient(data);
  GEVMLEHessianFunctor h;
  GEVMLEProblem mle_problem(mle, gradient, h);
  optimo::solvers::BFGS<double> bfgs;
  double min_val;
  auto res = bfgs(mle_problem, &x, &min_val);
  *mu = x(0);
  *sigma = x(1);
  *xi = x(2);
  exit_flag = res == 0;
  LOG_IF(INFO, !exit_flag) << "MLE did not converge: res=" << res;
  VLOG_IF(1, !exit_flag) << "Values: " << x.block<3, 1>(0, 0).transpose()
                         << " min_val: " << min_val;

  return exit_flag;
}
}  // evd
}  // distributions
}  // statx
