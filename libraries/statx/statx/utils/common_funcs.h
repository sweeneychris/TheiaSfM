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

#ifndef STATX_UTILS_COMMON_FUNCS_H_
#define STATX_UTILS_COMMON_FUNCS_H_

#include <Eigen/Core>
#include <cmath>
#include <vector>

namespace statx {
namespace utils {

using Eigen::Matrix;
using Eigen::Dynamic;
using std::vector;

/// Computes the sample mean
static inline
double mean(const vector<double>& y) {
  if (y.empty()) return 0.0;
  double mean = 0.0;
  for (double x : y) {
    mean += x;
  }
  return mean/y.size();
}

/// Computes the sample standard deviation
static inline
double stddev(const vector<double>& y) {
  if (y.empty()) return 0.0;
  double mu = mean(y);
  double s = 0.0;
  for (double x : y) {
    s += (x - mu)*(x - mu);
  }
  return sqrt(s / (y.size() - 1));
}

/// Computes the sample standard deviation provided the mean sample
static inline
double stddev(const vector<double>& y,
              const double mu) {
  if (y.empty()) return 0.0;
  double s = 0.0;
  for (double x : y) {
    s += (x - mu)*(x - mu);
  }
  return sqrt(s / (y.size() - 1));
}

/// Computes the sample mean vector
template <typename Scalar> inline
bool mean(const vector<Matrix<Scalar, Dynamic, 1> >& samples,
          Matrix<Scalar, Dynamic, 1>* mean) {
  if (!mean) return false;
  mean->setConstant(static_cast<Scalar>(0.0));
  const Scalar nsamples = static_cast<Scalar>(samples.size());
  for (const Matrix<Scalar, Dynamic, 1>& z : samples) *mean += z / nsamples;
  return true;
}

/// Computes the sample covariance matrix
template <typename Scalar> inline
bool cov(const vector<Matrix<Scalar, Dynamic, 1> >& samples,
         const Matrix<Scalar, Dynamic, 1>& mean,
         Matrix<Scalar, Dynamic, Dynamic>* cov_mat) {
  if (!cov_mat) return false;
  const int dim = mean.rows();
  cov_mat->resize(dim, dim);
  cov_mat->setConstant(static_cast<Scalar>(0.0));
  const Scalar norm_factor = static_cast<Scalar>(samples.size() - 1);
  for (const Matrix<Scalar, Dynamic, 1>& z : samples) {
    *cov_mat += (z - mean)*(z - mean).transpose() / norm_factor;
  }
  return true;
}

/// Computes the Mahalanobis distance
template <typename Scalar> inline
Scalar mahalanobis_distance(const Matrix<Scalar, Dynamic, 1>& x,
                            const Matrix<Scalar, Dynamic, 1>& mean_vec,
                            const Matrix<Scalar, Dynamic, Dynamic>& covar_mat) {
  return static_cast<Scalar>(sqrt(x.dot(covar_mat.inverse()*x)));
}
}  // utils
}  // statx
#endif  // STATX_UTILS_COMMON_FUNCS_H_
