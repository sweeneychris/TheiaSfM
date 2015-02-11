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

#include <glog/logging.h>
#include <optimo/solvers/bfgs.h>
#include <limits>
#include "statx/distributions/evd/gpd_mle.h"
#include "statx/utils/common_funcs.h"

namespace statx {
namespace distributions {
namespace evd {

using Eigen::Matrix;
using Eigen::Dynamic;
using std::vector;

// Implementation of optimo objects
// See pg 80 from Stuart Coles.
double GPDMLEObjective::operator()(const Matrix<double, Dynamic, 1>& x) const {
  double acc1 = 0.0;
  double acc2 = 0.0;
  const double& xi = x(0);
  const double& sigma = x(1);
  const size_t k = data_.size();
  const double lambda_inv = 1.0/lambda_;  // For penalties
  if (sigma <= 0.0) return std::numeric_limits<double>::max();
  for (double z : data_) {
    const double log_arg = 1.0 + xi*z/sigma;
    if (log_arg < 0.0) return std::numeric_limits<double>::max();
    acc1 += log(log_arg);
    // Penalty term for log_arg < 0
    acc2 += (log_arg <= lambda_inv) ? log(lambda_*log_arg) : 0.0;
  }
  const double penalty1 =
      (sigma <= lambda_inv) ? -alpha_*log(lambda_*sigma) : 0.0;
  return k*log(sigma) + (1.0 + 1.0/xi)*acc1 + penalty1 - beta_*acc2;
}

// Implementation of GPDMLEGradient
void GPDMLEGradientFunctor::operator()(const Matrix<double, Dynamic, 1>& x,
                                       Matrix<double, Dynamic, 1>* g) const {
  double acc1 = 0.0;
  double acc2 = 0.0;
  double acc3 = 0.0;
  double acc4 = 0.0;
  double acc5 = 0.0;
  const double& xi = x(0);
  const double& sigma = x(1);
  const double sigma_inv = 1.0/sigma;
  const size_t k = data_.size();
  const double lambda_inv = 1.0/lambda_;
  for (double z : data_) {
    const double sigma_xi_z = sigma + xi*z;
    const double sigma_xi_z_inv = 1.0 / sigma_xi_z;
    const double log_arg = 1.0 + xi*z*sigma_inv;
    acc1 += xi*z*sigma_xi_z_inv;
    acc2 += z*sigma_xi_z_inv;
    acc3 += log(log_arg);
    // Calculate terms for penalties
    bool log_arg_flag = log_arg <= lambda_inv;
    acc4 += log_arg_flag ? xi*z*sigma_xi_z_inv: 0.0;
    acc5 += log_arg_flag ? z*sigma_xi_z_inv : 0.0;
  }
  const double xi_term = (1.0 + 1.0/xi);
  (*g)(0) = xi_term*acc2 - acc3/(xi*xi);
  (*g)(1) = k*sigma_inv - sigma_inv*xi_term*acc1;
  // Terms from the penalties
  double penalty1 = sigma <= lambda_inv ? -alpha_*sigma_inv : 0.0;
  (*g)(0) -= beta_*acc5;
  (*g)(1) += penalty1 + beta_*sigma_inv*acc4;
}

// Implementation of GPDMLEHessian (dummy implementation)
void
GPDMLEHessianFunctor::operator()(const Matrix<double, Dynamic, 1>& x,
                                 Matrix<double, Dynamic, Dynamic>* h) const {
  
}

// GPD Problem for xi != 0
class GPDMLEProblem : public optimo::ProblemLS<double> {
 public:
  // Constructor
  GPDMLEProblem(const GPDMLEObjective& obj,
                const GPDMLEGradientFunctor& g,
                const GPDMLEHessianFunctor& h) :
      optimo::ProblemLS<double>(obj, g, h) { }

  // Destructor
  virtual ~GPDMLEProblem(void) { }
};

// Fit a gpd distribution using an MLE estimator
// When xi > 0.5, MLE might not converge
bool gpdfit_mle(const vector<double>& data,
                double* xi,
                double* sigma) {
  if (!xi || !sigma) return false;
  bool exit_flag = false;
  const int n = 2;  // 2 Prams to estimate
  Matrix<double, Dynamic, 1> x(n);  // Unkown vector
  // Solve using BFGS
  GPDMLEObjective mle(data);
  GPDMLEGradientFunctor gradient(data);
  GPDMLEHessianFunctor h;
  GPDMLEProblem mle_problem(mle, gradient, h);

  const double mean = statx::utils::mean(data);
  const double temp = statx::utils::stddev(data, mean);
  const double var = temp*temp;
  const double mean_sqrd= mean*mean;

  // TODO(vfragoso): There could be cases, where the initial points
  // map to an infeasible solution. Either, come up with better heuristics
  // or device a mechanism to be robust to that.
  // according to EVIR's package (Initial params)
  // xi0 <- -0.5 * (((xbar * xbar)/s2) - 1)
  // beta0 <- 0.5 * xbar * (((xbar * xbar)/s2) + 1)
  x(0) = 0.5 * (1.0 - mean_sqrd/var);  // xi 
  x(1) = 0.5 * mean * (mean_sqrd/var + 1.0);  // sigma
  VLOG(1) << "x0 = " << x.transpose();

  // // Checking the initial point
  // double min_z = std::numeric_limits<double>::min();
  // for (double z : data) if (min_z > z) min_z = z;
  // if (x(1)/min_z >= x(0)) {  // Not a valid starting point
  //   VLOG(1) << "Invalid starting point";
  //   if (x(1) > 0.0) {
  //     // xi = \sigma / min_i(z_i) + \delta
  //     x(1) = (x(0) / min_z) + 1e-3;
  //   }
  // }

  optimo::solvers::BFGS<double> bfgs;
  double min_val;
  auto res = bfgs(mle_problem, &x, &min_val);
  *xi = x(0);  // xi
  *sigma = x(1);  // sigma
  exit_flag = res == 0;
  LOG_IF(INFO, !exit_flag) << "MLE did not converge: res=" << res;
  VLOG_IF(1, !exit_flag) << "Values: " << x.transpose()
                         << " min_val: " << min_val;
  return exit_flag;
}
}  // evd
}  // distributions
}  // statx
