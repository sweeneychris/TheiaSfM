// Copyright (C) 2013  Victor Fragoso <vfragoso@cs.ucsb.edu>
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

#ifndef OPTIMO_CORE_NUMERICAL_GRADIENT_H_
#define OPTIMO_CORE_NUMERICAL_GRADIENT_H_

#include <Eigen/Core>
#include "optimo/core/objects.h"
#include "optimo/core/objects_ls.h"

namespace optimo {
/// Numerical Gradient Computation using the Secant Method

/// This class computes the gradient using the secant method.
/// \f[
/// \nabla f_0(\mathbf{x}) \approx
/// \frac{f_0(\mathbf{x}_r) - f_0(\mathbf{x}_l)}{2*h}
/// \f]
///
/// where
///
/// \f{eqnarray*}{
/// \mathbf{x}_r & = & \mathbf{x} + h\mathbb{1} \\
/// \mathbf{x}_l & = & \mathbf{x} - h\mathbb{1} \\
/// \f}
/// and \f$\mathbb{1}\f$ is a vector of ones.
///
/// This class is recommended for problems that fit in the stack. If this
/// is a concern, use the SecantGradientFunctorLS class.
template <typename Scalar, uint n>
class SecantGradientFunctor : public GradientFunctor<Scalar, n> {
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  /// Constructor
  SecantGradientFunctor(const ObjectiveFunctor<Scalar, n>& obj,
                        const Scalar h) :
      GradientFunctor<Scalar, n>(), objective_(obj), h_(h) { }

  // Destructor
  virtual ~SecantGradientFunctor(void) { }

  /// Compute gradient
  virtual void
  operator()(const Eigen::Matrix<Scalar, n, 1>& x,
             Eigen::Matrix<Scalar, n, 1>* g) const {
    Scalar den = 2*h_;
    Eigen::Matrix<Scalar, n, 1> xl, xr;
    for (unsigned int i = 0; i < n; i++) {
      xl = x;
      xr = x;
      xr(i, 0) += h_;
      xl(i, 0) -= h_;
      (*g)(i, 0) = (objective_(xr) - objective_(xl))/den;
    }
  }

 protected:
  const ObjectiveFunctor<Scalar, n>& objective_;
  const Scalar h_;
};
/// Numerical Gradient Computation using the Secant Method

/// This class computes the gradient using the secant method.
/// \f[
/// \nabla f_0(\mathbf{x}) \approx
/// \frac{f_0(\mathbf{x}_r) - f_0(\mathbf{x}_l)}{2*h}
/// \f]
///
/// where
///
/// \f{eqnarray*}{
/// \mathbf{x}_r & = & \mathbf{x} + h\mathbb{1} \\
/// \mathbf{x}_l & = & \mathbf{x} - h\mathbb{1} \\
/// \f}
/// and \f$\mathbb{1}\f$ is a vector of ones.
template <typename Scalar>
class SecantGradientFunctorLS : public GradientFunctorLS<Scalar> {
 public:
  /// Constructor
  SecantGradientFunctorLS(const ObjectiveFunctorLS<Scalar>& obj,
                          const Scalar h) :
      GradientFunctorLS<Scalar>(), objective_(obj), h_(h) { }

  virtual ~SecantGradientFunctorLS(void) { }

  /// Compute gradient
  virtual void
  operator()(const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& x,
             Eigen::Matrix<Scalar, Eigen::Dynamic, 1>* g) const {
    const int n = x.rows();
    Scalar den = 2*h_;
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> xl(n), xr(n);
    for (unsigned int i = 0; i < n; i++) {
      xl = x;
      xr = x;
      xr(i, 0) += h_;
      xl(i, 0) -= h_;
      (*g)(i, 0) = (objective_(xr) - objective_(xl))/den;
    }
  }

 protected:
  const ObjectiveFunctorLS<Scalar>& objective_;
  const Scalar h_;
};
}  // optimo

#endif  // OPTIMO_CORE_NUMERICAL_GRADIENT_H_
