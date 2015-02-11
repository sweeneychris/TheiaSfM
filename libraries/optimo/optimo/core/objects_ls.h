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

#ifndef OPTIMO_CORE_OBJECTS_LS_H_
#define OPTIMO_CORE_OBJECTS_LS_H_

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <vector>

// m: Number of Constraints
// n: Number of Variables
typedef unsigned int uint;
namespace optimo {
using Eigen::Dynamic;
using Eigen::Matrix;
/// Objective functor

/// This functor evaluates the objective function to optimize
/// This function assumes the form:
/// \f$ f_0: \mathbb{R}^n \rightarrow \mathbb{R} \f$
template <typename Scalar>
struct ObjectiveFunctorLS {
  /// Evaluates objective function
  virtual Scalar
  operator()(const Matrix<Scalar, Dynamic, 1>& x) const = 0;
};

/// Inequality constraints functor

/// This functor evaluates the constraint function \f$f_i(\mathbf{x})\f$
/// This function assumes the form:
/// \f$ f_i: \mathbb{R}^n \rightarrow \mathbb{R} \f$
template <typename Scalar>
struct ConstraintFunctorLS {
  /// Evaluates constraint
  virtual Scalar
  operator()(const Matrix<Scalar, Dynamic, 1>& x) const = 0;
};

/// Inequality constraints functor

/// This functor evaluates the constraint function \f$f_i(\mathbf{x})\f$
/// This function assumes the form:
/// \f$ f_i: \mathbb{R}^n \rightarrow \mathbb{R} \f$
template <typename Scalar>
struct ConstraintsLS {
  /// Returns a ConstraintFunctorLS reference
  virtual ConstraintFunctorLS<Scalar>&
  operator[](const unsigned int i) = 0;
};

/// Gradient functor

/// This functor evaluates the gradient of the objective function
/// \f$ \nabla f_0 (\mathbf{x})\f$
/// Assumes the form:
/// \f$ \nabla f_0 (\mathbf{x}): \mathbb{R}^n \rightarrow \mathbb{R}^n \f$
template <typename Scalar>
struct GradientFunctorLS {
  /// Computes the gradient
  virtual void
  operator()(const Matrix<Scalar, Dynamic, 1>& x,
             Matrix<Scalar, Dynamic, 1>* g) const = 0;
};

/// Dense Hessian functor

/// This functor evaluates the Hessian of the objective function
/// \f$ \nabla^2 f_0(\mathbf{x})\f$
/// Assumes the form:
/// \f[
/// \nabla^2 f_0(\mathbf{x}): \mathbb{R}^n \rightarrow \mathbb{R}^{n \times n}
/// \f]
template <typename Scalar>
struct HessianFunctorLS {
  /// Computes the Hessian
  virtual void 
  operator()(
      const Matrix<Scalar, Dynamic, 1>& x,
      Matrix<Scalar, Dynamic, Dynamic>* h) const = 0;
};

/// Dummy constraint functor for unconstrained problems.

/// Helper object for uncosntrained problems
template <typename Scalar>
struct VoidConstraintFunctorLS : public ConstraintFunctorLS<Scalar> {
  /// Dummy evaluation of a constraint functor. Returns 0 always.
  virtual Scalar
  operator()(const Matrix<Scalar, Dynamic, 1>& x) const {
    return static_cast<Scalar>(0.0);
  }
};

/// Helper object for unconstrained problem

/// This helper object is a singleton object.
template <typename Scalar>
class NoConstraintsLS : public ConstraintsLS<Scalar> {
 public:
  /// Returns the singleton instance.
  static NoConstraintsLS<Scalar>& getInstance(void);

  // Non copyable
  // Copy constructor
  NoConstraintsLS(const NoConstraintsLS<Scalar>& x) = delete;

  // Operator =
  NoConstraintsLS<Scalar>&
  operator=(const NoConstraintsLS<Scalar>& rhs) = delete;

  /// Dummy implementation for no constraints
  virtual ConstraintFunctorLS<Scalar>& operator[](const unsigned int i) {
    return *c;
  }

 protected:
  static VoidConstraintFunctorLS<Scalar>* c;

 private:
  static void cleanUp(void);

  // Constructor
  NoConstraintsLS(void) {
    c = new VoidConstraintFunctorLS<Scalar>;
    atexit(&cleanUp);
  }

  // Destructor
  ~NoConstraintsLS(void) { }

  static NoConstraintsLS<Scalar>* instance;
};

template <typename Scalar>
NoConstraintsLS<Scalar>* NoConstraintsLS<Scalar>::instance = nullptr;

template <typename Scalar>
VoidConstraintFunctorLS<Scalar>* NoConstraintsLS<Scalar>::c = nullptr;

// Cleanup function
template <typename Scalar>
void NoConstraintsLS<Scalar>::cleanUp(void) {
  delete c;
  c = nullptr;
}

template <typename Scalar>
NoConstraintsLS<Scalar>& NoConstraintsLS<Scalar>::getInstance(void) {
  if (!instance) {
    instance = new NoConstraintsLS<Scalar>;
  }
  return *instance;
}

/// Class defining a convex Problem

/// The convex problem is assumed to have the following form:
/// \f{eqnarray*}{
/// \underset{\mathbf{x}}{\text{minimize}} & f_0(\mathbf{x}) & \\
/// \text{subject to} & & \\
/// f_i(\mathbf{x}) & \leq & 0 \\
/// A_{\text{eq}} \mathbf{x} & = & \mathbf{b}_{\text{eq}}
/// \f}
///
/// This class Eigen matrices on the stack. Thus, use this class if
/// the problem fits in the stack.
template <typename Scalar>
class ProblemLS {
 public:
  const ObjectiveFunctorLS<Scalar>& objective;  ///< Objective functor
  const ConstraintsLS<Scalar>& constraints;  ///< Constraints
  const GradientFunctorLS<Scalar>& gradient;  ///< Gradient functor
  const HessianFunctorLS<Scalar>& hessian;  ///< Hessian functor
  const Matrix<Scalar, Dynamic, Dynamic> A;  //< Equality constraints matrix
  const Matrix<Scalar, Dynamic, 1> b;  ///< Equality vector

  /// Constrained Problem
  ProblemLS(const ObjectiveFunctorLS<Scalar>& obj,
          const ConstraintsLS<Scalar>& con,
          const Matrix<Scalar, Dynamic, Dynamic>& Aeq,
          const Matrix<Scalar, Dynamic, 1>& beq,
          const GradientFunctorLS<Scalar>& g,
          const HessianFunctorLS<Scalar>& h) :
      objective(obj), constraints(con), gradient(g), hessian(h), A(Aeq), b(beq)
  { }

  /// Unconstrained problem
  ProblemLS(const ObjectiveFunctorLS<Scalar>& obj,
            const GradientFunctorLS<Scalar>& g,
            const HessianFunctorLS<Scalar>& h) :
      objective(obj),
      constraints(NoConstraintsLS<Scalar>::getInstance()),
      gradient(g), hessian(h) { }

  // Destructor
  virtual ~ProblemLS(void) { }
};
}  // optimo
#endif  // OPTIMO_CORE_OBJECTS_LS_H_
