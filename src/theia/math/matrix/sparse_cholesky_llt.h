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

#ifndef THEIA_MATH_MATRIX_SPARSE_CHOLESKY_LLT_H_
#define THEIA_MATH_MATRIX_SPARSE_CHOLESKY_LLT_H_

#include <cholmod.h>
#include <Eigen/SparseCore>

// UF_long is deprecated but SuiteSparse_long is only available in
// newer versions of SuiteSparse. So for older versions of
// SuiteSparse, we define SuiteSparse_long to be the same as UF_long,
// which is what recent versions of SuiteSparse do anyways.
#ifndef SuiteSparse_long
#define SuiteSparse_long UF_long
#endif

namespace theia {

// A class for performing the choleksy decomposition of a sparse matrix using
// CHOLMOD from SuiteSparse. This allows us to utilize the supernodal algorithms
// which are not included with Eigen. CHOLMOD automatically determines if the
// simplicial or supernodal algorithm is the best choice. The interface is meant
// to mimic the Eigen linear solver interface except that it is not templated
// and requires sparse matrices.
//
// NOTE: The matrix mat should be a symmetric matrix.
class SparseCholeskyLLt {
 public:
  explicit SparseCholeskyLLt(const Eigen::SparseMatrix<double>& mat);
  SparseCholeskyLLt();
  ~SparseCholeskyLLt();

  // Perform symbolic analysis of the matrix. This is useful for analyzing
  // matrices with the same sparsity pattern when used in conjunction with
  // Factorize().
  void AnalyzePattern(const Eigen::SparseMatrix<double>& mat);

  // Perform numerical decomposition of the current matrix. If the matrix has
  // the same sparsity pattern as the previous decomposition then this method
  // may be used to efficiently decompose the matrix by avoiding symbolic
  // analysis.
  void Factorize(const Eigen::SparseMatrix<double>& mat);

  // Computes the Cholesky decomposition of mat. This is the same as calling
  // AnalyzePattern() followed by Factorize().
  void Compute(const Eigen::SparseMatrix<double>& mat);

  // Returns the current state of the decomposition. After each step users
  // should ensure that Info() returns Eigen::Success.
  Eigen::ComputationInfo Info();

  // Using the cholesky decomposition, solve for x that minimizes
  //    lhs * x = rhs
  // where lhs is the factorized matrix.
  Eigen::VectorXd Solve(const Eigen::VectorXd& rhs);

 private:
  cholmod_common cc_;
  cholmod_factor* cholmod_factor_;
  bool is_factorization_ok_, is_analysis_ok_;
  Eigen::ComputationInfo info_;
};

}  // namespace theia

#endif  // THEIA_MATH_MATRIX_SPARSE_CHOLESKY_LLT_H_
