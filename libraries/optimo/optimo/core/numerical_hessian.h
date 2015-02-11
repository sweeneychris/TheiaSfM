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

#ifndef OPTIMO_CORE_NUMERICAL_HESSIAN_H_
#define OPTIMO_CORE_NUMERICAL_HESSIAN_H_

#include <Eigen/Core>
#include "optimo/core/objects.h"

namespace optimo {
/// Numerical Hessian computation

/// This class calculates the Hessian matrix using the secant method.
template <typename Scalar, uint n>
class NumericalHessian : public HessianFunctor<Scalar, n> {
 public:
  /// Constructor
  NumericalHessian(const ObjectiveFunctor<Scalar, n>& obj,
                   const Scalar h) :
      HessianFunctor<Scalar, n>(), objective_(obj), h_(h) { }

  virtual ~NumericalHessian(void) { }

  /// Computes Hessian matrix
  virtual void
  operator()(const Eigen::Matrix<Scalar, n, 1>& x,
             Eigen::Matrix<Scalar, n, n>* H) const {
    Scalar den = h_*h_;
    // 1. Calculate Diagonal of the Hessian
    Eigen::Matrix<Scalar, n, 1> xnn, xpp, xnp, xpn;
    for (unsigned int i = 0; i < n; i++) {
      xnn = x;
      xpp = x;
      xpp(i, 0) += h_;
      xnn(i, 0) -= h_;
      (*H)(i, i) = (objective_(xpp) - 2*objective_(x) + objective_(xnn))/den;
    }

    // 2. Compute the upper triangular part of the Hessian
    den *= 4;
    for (unsigned int i = 0; i < n; i++) {
      for (unsigned int j = i + 1; j < n; j++) {
        xpp = x;
        xnp = x;
        xpn = x;
        xnn = x;

        xpp(i, 0) += h_;
        xpp(j, 0) += h_;

        xpn(i, 0) += h_;
        xpn(j, 0) -= h_;

        xnp(i, 0) -= h_;
        xnp(j, 0) += h_;

        xnn(i, 0) -= h_;
        xnn(j, 0) -= h_;

        (*H)(i, j) = (objective_(xpp) - objective_(xpn) -
                      objective_(xnp) + objective_(xnn))/den;
        (*H)(j, i) = (*H)(i, j);
      }
    }
  }

 protected:
  const ObjectiveFunctor<Scalar, n>& objective_;
  const Scalar h_;
};
}  // optimo
#endif  // OPTIMO_CORE_NUMERICAL_HESSIAN_H_
