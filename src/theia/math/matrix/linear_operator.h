// Copyright (C) 2017 The Regents of the University of California (Regents).
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
#include <glog/logging.h>

namespace theia {

// This is an abstract base class for linear operators. It supports
// access to size information and left and right multiply operators. This class
// is inspired by the LinearOperator class in the Ceres Solver library:
// www.ceres-solver.org
class LinearOperator {
 public:
  virtual ~LinearOperator()=default;

  // y = y + Ax;
  virtual void RightMultiply(const Eigen::VectorXd& x,
                             Eigen::VectorXd* y) const = 0;
  // y = y + A'x;
  virtual void LeftMultiply(const Eigen::VectorXd& x,
                            Eigen::VectorXd* y) const = 0;

  virtual int num_rows() const = 0;
  virtual int num_cols() const = 0;
};

}  // namespace theia

#endif  // THEIA_MATH_MATRIX_LINEAR_OPERATOR_H_
