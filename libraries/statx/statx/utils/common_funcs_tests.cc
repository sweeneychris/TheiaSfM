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

#include "statx/utils/common_funcs.h"
#include <glog/logging.h>
#include <random>
#include <cmath>
#include "gtest/gtest.h"

namespace statx {
namespace utils {
using Eigen::VectorXd;
using Eigen::MatrixXd;

// Sample generated with mvnrnd(zeros(4, 1), eye(4), 25)
const vector<double> mvn_sample {
  8.017041e-01, 3.830239e-01, 5.406331e-01, -1.281281e+00,
      1.053305e+00, 4.120354e-01, 9.758409e-01, -2.203264e+00,
      -7.488768e-01, 4.054926e-01, -1.568704e-01, -5.712463e-01,
      -9.363265e-01, -3.637808e-01, 2.777993e-01, 2.139965e-01,
      -1.269087e+00, -5.992720e-01, 6.395173e-01, 9.423769e-01,
      4.979806e-01, -5.895890e-01, -8.097802e-02, 9.372549e-02,
      2.789081e+00, 8.535408e-01, 5.408701e-01, -1.122312e+00,
      7.275720e-01, -1.853008e+00, -1.262565e+00, 3.061578e-01,
      -7.730641e-01, -2.073032e-01, 1.110424e+00, -1.172335e+00,
      8.366338e-01, 2.703782e-01, -9.895627e-01, -9.609666e-01,
      -1.128330e+00, -6.527710e-01, -1.828836e+00, -6.537350e-01,
      -1.424470e+00, 4.772273e-01, 1.384499e+00, -1.229394e+00,
      7.174423e-01, -7.131965e-02, -6.272679e-02, -2.709651e-01,
      -7.779055e-01, -9.383013e-01, 4.489211e-01, -8.999501e-01,
      3.159859e-01, 1.613635e-01, -3.632585e-01, -2.856861e-01,
      1.406535e+00, -2.681829e-01, -1.020583e+00, -4.624215e-01,
      4.011246e-01, -4.098726e-01, -3.072989e+00, -4.097852e-01,
      9.296603e-01, -7.113227e-01, 6.262790e-01, -5.035390e-01,
      -1.605802e+00, 6.144548e-02, -2.866845e-01, 1.233297e+00,
      6.615362e-01, -1.846129e+00, -1.973429e-01, 6.103052e-01,
      2.138502e+00, -3.983331e-01, 4.056054e-01, 5.907216e-02,
      5.411394e-01, -5.435481e-01, -1.419348e+00, -1.466947e+00,
      -1.540877e+00, -9.118985e-01, -7.294452e-01, -1.625803e+00,
      -2.031428e-01, 6.526986e-01, 1.147328e+00, -1.964752e+00,
      -4.999652e-01, -7.342713e-01, 5.978646e-01, 2.605196e+00 };

TEST(CommonFuncs, Mean) {
  vector<double> samples {2, 2, 2, 2, 2, 2, 2};
  const double m = mean(samples);
  const double m_gt = 2.0;
  ASSERT_NEAR(m_gt, m, 1e-3);
}

TEST(CommonFuncs, Stddev) {
  vector<double> samples {2.1, 2.2, 1.8, 1.9, 2.0};
  const double s = stddev(samples);
  const double s_gt = 0.1581;
  ASSERT_NEAR(s_gt, s, 1e-3);
  const double mu = mean(samples);
  const double s2 = stddev(samples, mu);
  ASSERT_NEAR(s_gt, s, 1e-3);
}

TEST(CommonFuncs, VecMean) {
  const int dim = 4;
  const int nsamples = 25;
  VectorXd mean_vec(dim);
  VectorXd mean_gt = VectorXd::Zero(dim);
  vector<VectorXd> vec_samples(nsamples, VectorXd::Zero(dim));
  int k = 0;
  for (int i = 0; i < nsamples; i++) {
    for (int j = 0; j < dim; j++) {
      vec_samples[i][j] = mvn_sample[k++];
    }
  }
  bool exit_flag = mean(vec_samples, &mean_vec);
  VLOG(1) << mean_vec.transpose();
  ASSERT_TRUE(exit_flag);
  double error = (mean_gt - mean_vec).norm();
  ASSERT_LT(error, 1.0);
}

TEST(CommonFuncs, CovarianceEstimation) {
  const int dim = 4;
  const int nsamples = 25;
  MatrixXd cov_mat(dim, dim);
  MatrixXd cov_mat_gt(dim, dim);
  cov_mat_gt << 1.3495, 0.0952, -0.0745, -0.2005,
      0.0952, 0.4737, 0.2495, -0.3260,
      -0.0745, 0.2495, 1.1056, -0.0800,
      -0.2005, -0.3260, -0.0800, 1.1432;
  VectorXd mean_vec(dim);
  VectorXd mean_gt = VectorXd::Zero(dim);
  vector<VectorXd> vec_samples(nsamples, VectorXd::Zero(dim));
  int k = 0;
  for (int i = 0; i < nsamples; i++) {
    for (int j = 0; j < dim; j++) {
      vec_samples[i][j] = mvn_sample[k++];
    }
  }
  mean(vec_samples, &mean_vec);
  ASSERT_TRUE(cov(vec_samples, mean_vec, &cov_mat));
  const double error = (cov_mat - cov_mat_gt).norm();
  ASSERT_LT(error, 1.0);
}
}  // utils
}  // statx
