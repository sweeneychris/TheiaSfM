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

#include "statx/distributions/evd/gev_ceres.h"
#include "statx/utils/common_funcs.h"
#include "statx/utils/ecdf.h"

namespace statx {
namespace distributions {
namespace evd {

using ceres::CostFunction;
using ceres::Problem;
using ceres::Solver;
using ceres::Solve;
using std::vector;

bool GEVCostFunctionAnalytic::Evaluate(double const* const* parameters,
                                       double* residuals,
                                       double** jacobians) const {
  const double& mu = parameters[0][0];  // location
  const double& sigma = parameters[0][1];  // scale
  const double& xi = parameters[0][2];  // shape
  const double log_term = -log(p_);
  const double pow_term = pow(log_term, -xi);
  const double xi_term = (1.0 - pow_term)/xi;

  // Residual function
  // f(\theta) = z - \mu + (\sigma/\xi)*(1 - (-\log(p))^{-\xi})
  // where \theta = [\mu \sigma \xi]'
  residuals[0] = z_ - mu + sigma*xi_term;

  if (jacobians && jacobians[0]) {
    // Jacobian of f(\theta)
    // df/d\mu = -1
    // df/d\sigma = (\frac{1 - (-\log(p))^{-\xi}}{\xi})
    // df/d\xi = term1 + term2
    // term1 = -(1 - (-\log(p))^{-\xi})*sigma/\xi^{-2}
    // term2 = (\sigma/\xi)*(-\log(p))^{-\xi} \log((-\log(p))^{-\xi})
    const double sigma_xi = sigma/xi;
    const double term1 = -xi_term*sigma_xi;
    const double term2 = sigma_xi*pow_term*log(log_term);
    jacobians[0][0] = -1;
    jacobians[0][1] = xi_term;
    jacobians[0][2] = term1 + term2;
  }
  return true;
}

// GEV fit using the quantile least-squares method
bool gevfit_ceres(const vector<double>& data,
                  double* mu,
                  double* sigma,
                  double* xi) {
  // TODO(vfragoso): Should I consider avoiding temporal variables?
  const double var = statx::utils::stddev(data);
  const double sigma0 = sqrt(6*var*var)/M_PI;  // sigma
  const double mu0 = statx::utils::mean(data) -0.57722*sigma0;  // mu
  const double xi0 = 0.1;  // according to EVIR's package

  // Initialize Params
  *mu = mu0;
  *sigma = sigma0;
  *xi = xi0;

  VLOG(1) << "Initial params: mu=" << mu0 << " sigma=" << sigma0
          << " xi=" << xi0;
  // Calculate the ECDF from the data
  vector<double> fx, x;
  statx::utils::ecdf(data, &fx, &x);

  // Build Ceres Objects
  Problem problem;
  // Only to the second to last element, as the last element has the
  // percentile of 1, and that causes numerical problems.
  for (int i = 0; i < fx.size() - 1; i++) {
    CostFunction* cost = new GEVCostFunctionAnalytic(x[i], fx[i]);
    problem.AddResidualBlock(cost, NULL, mu, sigma, xi);
  }

  // Solve!
  Solver::Options options;
  options.max_num_iterations = 100;
  options.linear_solver_type = ceres::DENSE_QR;
  options.minimizer_progress_to_stdout = false;
  options.update_state_every_iteration = true;

  // Get Summary!
  Solver::Summary summary;
  Solve(options, &problem, &summary);
  VLOG(3) << "\n" << summary.BriefReport();
  auto term_type = summary.termination_type;
  bool exit_flag =
      term_type == ceres::CONVERGENCE ||
      term_type == ceres::USER_SUCCESS;

  return exit_flag;
}
}  // evd
}  // distributions
}  // statx
