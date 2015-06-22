// Copyright (C) 2015 The Regents of the University of California (Regents).
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

#ifndef THEIA_MATH_MATRIX_LINEAR_OPERATOR_H_
#define THEIA_MATH_MATRIX_LINEAR_OPERATOR_H_

#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <glog/logging.h>

namespace theia {

// A pure virtual class that will specify multiply methods. This will allow
// custom implementation of multiplication e.g., using sparse matrices or linear
// solving.
//
// NOTE: This class was inspired from the LinearOperator class of Ceres Solver:
// http://www.ceres-solver.org
class LinearOperator {
 public:
  virtual ~LinearOperator() {}

  // y = A*x
  virtual void RightMultiply(const Eigen::VectorXd& x,
                             Eigen::VectorXd* y) const = 0;

  virtual int Cols() const = 0;
  virtual int Rows() const = 0;
};

// A standard linear operator for dense matrices.
class DenseLinearOperator : public LinearOperator{
 public:
  explicit DenseLinearOperator(const Eigen::MatrixXd& A) : A_(A) {}

  // y = A*x
  virtual void RightMultiply(const Eigen::VectorXd& x,
                             Eigen::VectorXd* y) const {
    *y = A_ * x;
  }

  virtual int Cols() const { return A_.cols(); }
  virtual int Rows() const { return A_.rows(); }

 private:
  const Eigen::MatrixXd& A_;
};

// A standard linear operator for sparse matrices.
class SparseLinearOperator : public LinearOperator{
 public:
  explicit SparseLinearOperator(const Eigen::SparseMatrix<double>& A) : A_(A) {}

  // y = A*x
  virtual void RightMultiply(const Eigen::VectorXd& x,
                             Eigen::VectorXd* y) const {
    *y = A_ * x;
  }

  virtual int Cols() const { return A_.cols(); }
  virtual int Rows() const { return A_.rows(); }

 private:
  const Eigen::SparseMatrix<double>& A_;
};

// An inverse linear operator that can be used with power iterations to
// determine the smallest eigenvalues and eigenvectors of a matrix A.
class DenseInverseLULinearOperator : public LinearOperator {
 public:
  explicit DenseInverseLULinearOperator(const Eigen::MatrixXd& A) : A_(A) {}

  virtual void RightMultiply(const Eigen::VectorXd& x,
                             Eigen::VectorXd* y) const {
    *y = A_.fullPivLu().solve(x);
  }

  virtual int Cols() const { return A_.cols(); }
  virtual int Rows() const { return A_.rows(); }

 private:
  const Eigen::MatrixXd& A_;
};

// An inverse linear operator that can be used with power iterations to
// determine the smallest eigenvalues and eigenvectors of a matrix A.
class SparseInverseLULinearOperator : public LinearOperator{
 public:
  explicit SparseInverseLULinearOperator(const Eigen::SparseMatrix<double>& A)
      : A_(A) {
    linear_solver_.compute(A);
    CHECK_EQ(linear_solver_.info(), Eigen::Success)
        << "Sparse LU Decomposition failed.";
  }

  virtual void RightMultiply(const Eigen::VectorXd& x,
                             Eigen::VectorXd* y) const {
    *y = linear_solver_.solve(x);
  }

  virtual int Cols() const { return A_.cols(); }
  virtual int Rows() const { return A_.rows(); }

 private:
  const Eigen::SparseMatrix<double>& A_;
  Eigen::SparseLU<Eigen::SparseMatrix<double> > linear_solver_;
};

}  // namespace theia

#endif  // THEIA_MATH_MATRIX_LINEAR_OPERATOR_H_
