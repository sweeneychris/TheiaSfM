// Copyright (C) 2016 The Regents of the University of California (Regents).
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
//     * Neither the name of The Regents or University of California nor the
//       names of its contributors may be used to endorse or promote products
//       derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Please contact the author of this library if you have any questions.
// Author: Chris Sweeney (cmsweeney@cs.ucsb.edu)

// **  NOTE: Parts of this file were borrowed or inspired from the Ceres   **
// ** Solver library. The Ceres license is included here for completeness. **
//
//
// Ceres Solver - A fast non-linear least squares minimizer
// Copyright 2015 Google Inc. All rights reserved.
// http://ceres-solver.org/
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
// * Neither the name of Google Inc. nor the names of its contributors may be
//   used to endorse or promote products derived from this software without
//   specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Author: sameeragarwal@google.com (Sameer Agarwal)

#include "theia/math/matrix/sparse_cholesky_llt.h"

#include <cholmod.h>
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <glog/logging.h>

// UF_long is deprecated but SuiteSparse_long is only available in
// newer versions of SuiteSparse. So for older versions of
// SuiteSparse, we define SuiteSparse_long to be the same as UF_long,
// which is what recent versions of SuiteSparse do anyways.
#ifndef SuiteSparse_long
#define SuiteSparse_long UF_long
#endif

namespace theia {
namespace {

cholmod_sparse ViewAsCholmod(const Eigen::SparseMatrix<double>& const_mat) {
  Eigen::SparseMatrix<double>& mat =
      const_cast<Eigen::SparseMatrix<double>&>(const_mat);
  cholmod_sparse res;
  res.nzmax = mat.nonZeros();
  res.nrow = mat.rows();
  res.ncol = mat.cols();
  res.p = mat.outerIndexPtr();
  res.i = mat.innerIndexPtr();
  res.x = mat.valuePtr();
  res.sorted = 1;
  if (mat.isCompressed()) {
    res.packed = 1;
  } else {
    res.packed = 0;
    res.nz = mat.innerNonZeroPtr();
  }

  // Set to 0 if the matrix is not symmetric.
  res.stype = 1;
  res.itype = CHOLMOD_INT;
  res.xtype = CHOLMOD_REAL;
  res.dtype = CHOLMOD_DOUBLE;
  return res;
}

cholmod_dense ViewAsCholmod(const Eigen::VectorXd& const_vec) {
  Eigen::VectorXd& vec = const_cast<Eigen::VectorXd&>(const_vec);
  cholmod_dense res;
  res.nrow = vec.rows();
  res.ncol = vec.cols();
  res.nzmax = res.nrow * res.ncol;
  res.d = vec.size();
  res.x = reinterpret_cast<void*>(vec.data());
  res.z = 0;
  res.xtype = CHOLMOD_REAL;
  return res;
}
}  // namespace

// A class for performing the choleksy decomposition of a sparse matrix using
// CHOLMOD from SuiteSparse. This allows us to utilize the supernodal algorithms
// which are not included with Eigen. The interface is meant to mimic the Eigen
// linear solver interface except that it is not templated and requires sparse
// matrices.
SparseCholeskyLLt::SparseCholeskyLLt(
    const Eigen::SparseMatrix<double>& mat)
    : cholmod_factor_(nullptr),
      is_factorization_ok_(false),
      is_analysis_ok_(false),
      info_(Eigen::Success) {
  cholmod_start(&cc_);
  Compute(mat);
}

SparseCholeskyLLt::SparseCholeskyLLt()
    : cholmod_factor_(nullptr),
      is_factorization_ok_(false),
      is_analysis_ok_(false),
      info_(Eigen::Success) {
  cholmod_start(&cc_);
}

SparseCholeskyLLt::~SparseCholeskyLLt() {
  if (cholmod_factor_ != nullptr) {
    cholmod_free_factor(&cholmod_factor_, &cc_);
  }
  cholmod_finish(&cc_);
}

void SparseCholeskyLLt::AnalyzePattern(
    const Eigen::SparseMatrix<double>& mat) {
  // Release the current decomposition if there is one.
  if (cholmod_factor_ != nullptr) {
    cholmod_free_factor(&cholmod_factor_, &cc_);
  }

  // Get the cholmod view of the sparse matrix.
  cholmod_sparse A = ViewAsCholmod(mat);

  // Cholmod can try multiple re-ordering strategies to find a fill
  // reducing ordering. Here we just tell it use AMD with automatic
  // matrix dependence choice of supernodal versus simplicial
  // factorization.
  cc_.nmethods = 1;
  cc_.method[0].ordering = CHOLMOD_AMD;
  cc_.supernodal = CHOLMOD_AUTO;

  // Perform symbolic analysis of the matrix.
  cholmod_factor_ = cholmod_analyze(&A, &cc_);
  if (VLOG_IS_ON(2)) {
    cholmod_print_common(const_cast<char*>("Symbolic Analysis"), &cc_);
  }

  if (cc_.status != CHOLMOD_OK) {
    VLOG(2) << "cholmod_analyze failed. error code: %d" << cc_.status;
    info_ = Eigen::NumericalIssue;
    return;
  }

  is_analysis_ok_ = true;
  is_factorization_ok_ = false;
  info_ = Eigen::Success;
}

void SparseCholeskyLLt::Factorize(const Eigen::SparseMatrix<double>& mat) {
  CHECK_NOTNULL(cholmod_factor_);
  CHECK(is_analysis_ok_) << "Cannot call Factorize() because symbolic analysis "
                            "of the matrix (i.e. AnalyzePattern()) failed!";
  cholmod_sparse A = ViewAsCholmod(mat);

  // Save the current print level and silence CHOLMOD, otherwise
  // CHOLMOD is prone to dumping stuff to stderr, which can be
  // distracting when the error (matrix is indefinite) is not a fatal
  // failure.
  const int old_print_level = cc_.print;
  cc_.print = 0;

  cc_.quick_return_if_not_posdef = 1;
  int cholmod_status = cholmod_factorize(&A, cholmod_factor_, &cc_);
  cc_.print = old_print_level;

  // TODO(sameeragarwal): This switch statement is not consistent. It
  // treats all kinds of CHOLMOD failures as warnings. Some of these
  // like out of memory are definitely not warnings. The problem is
  // that the return value Cholesky is two valued, but the state of
  // the linear solver is really three valued. SUCCESS,
  // NON_FATAL_FAILURE (e.g., indefinite matrix) and FATAL_FAILURE
  // (e.g. out of memory).
  switch (cc_.status) {
    case CHOLMOD_NOT_INSTALLED:
      LOG(ERROR) << "CHOLMOD failure: Method not installed.";
      info_ = Eigen::InvalidInput;
      return;
    case CHOLMOD_OUT_OF_MEMORY:
      LOG(ERROR) << "CHOLMOD failure: Out of memory.";
      info_ = Eigen::NumericalIssue;
      return;
    case CHOLMOD_TOO_LARGE:
      LOG(ERROR) << "CHOLMOD failure: Integer overflow occured.";
      info_ = Eigen::NumericalIssue;
      return;
    case CHOLMOD_INVALID:
      LOG(ERROR) << "CHOLMOD failure: Invalid input.";
      info_ = Eigen::InvalidInput;
      return;
    case CHOLMOD_NOT_POSDEF:
      LOG(ERROR) << "CHOLMOD warning: Matrix not positive definite.";
      info_ = Eigen::NumericalIssue;
      return;
    case CHOLMOD_DSMALL:
      LOG(ERROR) << "CHOLMOD warning: D for LDL' or diag(L) or "
                    "LL' has tiny absolute value.";
      info_ = Eigen::NumericalIssue;
      return;
    case CHOLMOD_OK:
      // If everything is done successfully, set the appropriate flags to
      // success and exit.
      if (cholmod_status != 0) {
        info_ = Eigen::Success;
        is_factorization_ok_ = true;
        return;
      }
      LOG(ERROR) << "CHOLMOD failure: cholmod_factorize returned false "
                    "but cholmod_common::status is CHOLMOD_OK.";
      info_ = Eigen::NumericalIssue;
      return;
    default:
      LOG(ERROR) << "Unknown cholmod return code: " << cc_.status;
      info_ = Eigen::InvalidInput;
      return;
  }
}

void SparseCholeskyLLt::Compute(const Eigen::SparseMatrix<double>& mat) {
  AnalyzePattern(mat);
  Factorize(mat);
}

Eigen::ComputationInfo SparseCholeskyLLt::Info() {
  return info_;
}

// Using the cholesky decomposition, solve for x that minimizes
//    lhs * x = rhs
// where lhs is the factorized matrix.
Eigen::VectorXd SparseCholeskyLLt::Solve(const Eigen::VectorXd& rhs) {
  CHECK_NOTNULL(cholmod_factor_);
  CHECK(is_analysis_ok_) << "Cannot call Solve() because symbolic analysis "
                            "of the matrix (i.e. AnalyzePattern()) failed!";
  CHECK(is_factorization_ok_)
      << "Cannot call Solve() because numeric factorization "
         "of the matrix (i.e. Factorize()) failed!";

  Eigen::VectorXd solution;
  if (cc_.status != CHOLMOD_OK) {
    LOG(ERROR) << "cholmod_solve failed. CHOLMOD status is not CHOLMOD_OK";
    return solution;
  }

  // returns a cholmod_dense*
  cholmod_dense b = ViewAsCholmod(rhs);
  cholmod_dense* x = cholmod_solve(CHOLMOD_A, cholmod_factor_, &b, &cc_);
  if (x == nullptr) {
    info_ = Eigen::NumericalIssue;
    return solution;
  }
  solution = Eigen::Map<Eigen::VectorXd>(reinterpret_cast<double*>(x->x),
                                         x->nrow,
                                         x->ncol);
  return solution;
}

}  // namespace theia
