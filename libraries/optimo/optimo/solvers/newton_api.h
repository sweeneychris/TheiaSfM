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

#ifndef OPTIMO_SOLVERS_NEWTON_API_H_
#define OPTIMO_SOLVERS_NEWTON_API_H_
#include "optimo/core/objects.h"
#include "optimo/utils/matrix_utils.h"
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/Cholesky>
#include "optimo/solvers/solver.h"

namespace optimo {
namespace solvers {
/// Implements Newton's algorithm

/// Implementation of a general Newton's algorithm. This algorithm can solve
/// unconstrained as well as equality constrained problem. This implementation
/// is recommended for small size Problems that fit in the stack. 
/// For more info, please consult Convex Optimization by Stephen Boyd and
/// Lieven Vandenberghe.
template <typename Scalar, uint m, uint n>
class Newton : public Solver<Scalar> {
 public:
  enum Type {
    UNCONSTRAINED,
    EQUALITY_CONSTRAINED,
    INFEASIBLE_EQUALITY_CONSTRAINED
  };

  /// Constructor for Newton Solver
  Newton(const Type type = UNCONSTRAINED /**< Problem type */) : type_(type) {
    v_.setConstant(static_cast<Scalar>(1.0));
    switch (type) {
      case UNCONSTRAINED:
        method_ = new UnconstrainedMethod(this->options.alpha_,
                                          this->options.beta_);
        break;
      case EQUALITY_CONSTRAINED:
        method_ = new EqualityConstrainedMethod(this->options.alpha_,
                                                this->options.beta_);
        break;
      case INFEASIBLE_EQUALITY_CONSTRAINED:
        method_ = new InfeasibleMethod(this->options.alpha_,
                                       this->options.beta_, &v_);
        break;
    }
  }

  // Destructor
  virtual ~Newton(void) {
    delete method_;
  }

  /// Solve the problem
  TERMINATION_TYPE
  operator()(const Problem<Scalar, m, n>& problem,
             Eigen::Matrix<Scalar, n, 1>* x,
             Scalar* min_value);

 protected:
  //////////////////////////
  // Type (Method)
  struct Method {
    virtual ~Method(void) { }

    virtual TERMINATION_TYPE
    operator()(const Problem<Scalar, m, n>& problem,
               const Eigen::Matrix<Scalar, n, n>& hessian,
               const Eigen::Matrix<Scalar, n, 1>& gradient,
               const Scalar& epsilon,
               Eigen::Matrix<Scalar, n, 1>* x,
               Scalar* min_value) = 0;

    const Scalar alpha_;
    const Scalar beta_;

    Method(const Scalar alpha, const Scalar beta) :
        alpha_(alpha), beta_(beta) { }
    Eigen::Matrix<Scalar, n, 1> xnt_;  // Newton Step (Dense)
  };

  struct UnconstrainedMethod : public Method {
    UnconstrainedMethod(const Scalar alpha, const Scalar beta) :
        Method(alpha, beta) { }

    virtual TERMINATION_TYPE
    operator()(const Problem<Scalar, m, n>& problem,
               const Eigen::Matrix<Scalar, n, n>& hessian,
               const Eigen::Matrix<Scalar, n, 1>& gradient,
               const Scalar& epsilon,
               Eigen::Matrix<Scalar, n, 1>* x,
               Scalar* min_value);

    virtual
    Scalar line_search(const ObjectiveFunctor<Scalar, n>& objective,
                       const Eigen::Matrix<Scalar, n, 1>& x,
                       const Eigen::Matrix<Scalar, n, 1>& g);

    Scalar t_;
    Scalar fx_;
    Scalar lambda_sqrd_;
  };

  struct EqualityConstrainedMethod : public UnconstrainedMethod {
    EqualityConstrainedMethod(const Scalar alpha, const Scalar beta) :
        UnconstrainedMethod(alpha, beta) { }

    // Overriding
    virtual TERMINATION_TYPE
    operator()(const Problem<Scalar, m, n>& problem,
               const Eigen::Matrix<Scalar, n, n>& hessian,
               const Eigen::Matrix<Scalar, n, 1>& gradient,
               const Scalar& epsilon,
               Eigen::Matrix<Scalar, n, 1>* x,
               Scalar* min_value);
  };

  struct InfeasibleMethod : public Method {
    InfeasibleMethod(const Scalar alpha,
                     const Scalar beta,
                     Eigen::Matrix<Scalar, m, 1>* v) :
        Method(alpha, beta), v_(v) { }

    virtual TERMINATION_TYPE
    operator()(const Problem<Scalar, m, n>& problem,
               const Eigen::Matrix<Scalar, n, n>& hessian,
               const Eigen::Matrix<Scalar, n, 1>& gradient,
               const Scalar& epsilon,
               Eigen::Matrix<Scalar, n, 1>* x,
               Scalar* min_value);

    Scalar line_search(const ObjectiveFunctor<Scalar, n>& objective,
                       const Eigen::Matrix<Scalar, n, 1>& x,
                       const Eigen::Matrix<Scalar, m, 1>& v,
                       const Eigen::Matrix<Scalar, m, n>& A,
                       const Eigen::Matrix<Scalar, m, 1>& b,
                       const Eigen::Matrix<Scalar, n, 1>& g);

    Scalar t_;
    Eigen::Matrix<Scalar, m, 1>* v_;  // Dual Eq.
    Eigen::Matrix<Scalar, m, 1> vnt_;  // Dual Eq. step # of rows is m
  };

  /////////////////////////
  // Members
  // const Scalar epsilon;
  const Type type_;
  // const int max_iter;
  Method* method_;
  Eigen::Matrix<Scalar, m, 1> v_;  // Ineq. Constraints Lagrange mult.
};
}  // solvers
}  // optimo
#endif  // OPTIMO_SOLVERS_NEWTON_API_H_
