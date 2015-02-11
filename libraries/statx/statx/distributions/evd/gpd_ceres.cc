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

#include "statx/distributions/evd/gpd_ceres.h"
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

bool GPDCostFunctionAnalytic::Evaluate(double const* const* parameters,
                                       double* residuals,
                                       double** jacobians) const {
  const double& sigma = parameters[0][0];  // location
  const double& xi = parameters[0][1];  // scale

  // Resdual function
  // f(\theta) = z - sigma*((1 - p)^-xi - 1)/xi
  const double k = pow(1 - p_, -xi) - 1.0;
  residuals[0] = z_ - sigma * k / xi;
  
  if (jacobians && jacobians[0]) {
    // Jacobian of f(\theta)
    // df/d\sigma = -k / \sigma
    // df/d\xi = \sigma * k / \xi^2
    jacobians[0][0] = -k / xi;
    jacobians[0][1] = -jacobians[0][0] * sigma / xi;
  }
  return true;
}

bool gpdfit_ceres(const vector<double>& data,
                  double* sigma,
                  double* xi) {
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

  // TODO(vfragoso): Should I consider avoiding temporal variables?
  const double xi0 = 0.5 * (1.0 - mean_sqrd/var);  // xi 
  const double sigma0 = 0.5 * mean * (mean_sqrd/var + 1.0);  // sigma

  // Initialize Params
  *xi = xi0;
  *sigma = sigma0;

  VLOG(1) << "Initial params: sigma0=" << sigma0 << " xi0=" << xi0;

  // Calculate the ECDF from the data
  vector<double> fx, x;
  statx::utils::ecdf(data, &fx, &x);

  // Build Ceres Objects
  Problem problem;
  // Only to the second to last element, as the last element has the
  // percentile of 1, and that causes numerical problems.
  for (int i = 0; i < fx.size() - 1; i++) {
    CostFunction* cost = new GPDCostFunctionAnalytic(x[i], fx[i]);
    problem.AddResidualBlock(cost, NULL, sigma, xi);
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
