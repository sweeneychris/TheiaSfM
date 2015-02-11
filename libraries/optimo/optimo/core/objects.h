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

#ifndef OPTIMO_CORE_OBJECTS_H_
#define OPTIMO_CORE_OBJECTS_H_

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <vector>

#define NO_CONSTRAINTS 1

// m: Number of Constraints
// n: Number of Variables
typedef unsigned int uint;
namespace optimo {
/////////////////////////////////////////////////
// Small problems that fit in the Stack (n<16, m<16)

/// Objective functor

/// This functor evaluates the objective function to optimize
/// This function assumes the form:
/// \f$ f_0: \mathbb{R}^n \rightarrow \mathbb{R} \f$
template <typename Scalar, uint n>
struct ObjectiveFunctor {
  /// Evaluates objective function
  virtual Scalar operator()(const Eigen::Matrix<Scalar, n, 1>& x) const = 0;
};

/// Inequality constraints functor

/// This functor evaluates the constraint function \f$f_i(\mathbf{x})\f$
/// This function assumes the form:
/// \f$ f_i: \mathbb{R}^n \rightarrow \mathbb{R} \f$
template <typename Scalar, uint n>
struct ConstraintFunctor {
  /// Evaluates constraint function
  virtual Scalar operator()(const Eigen::Matrix<Scalar, n, 1>& x) const = 0;
};

/// Constraints object

/// This object collects all the constraint functions considered in the
/// optimization problem. The object assumes we have \f$m\f$ constraints.
template <typename Scalar, uint m, uint n>
struct Constraints {
  /// Returns a ConstraintFunctor
  virtual ConstraintFunctor<Scalar, n>&
  operator[](const unsigned int i) = 0;
};

/// Gradient functor

/// This functor evaluates the gradient of the objective function
/// \f$ \nabla f_0 (\mathbf{x})\f$
/// Assumes the form:
/// \f$ \nabla f_0 (\mathbf{x}): \mathbb{R}^n \rightarrow \mathbb{R}^n \f$
template <typename Scalar, uint n>
struct GradientFunctor {
  /// Computes the gradient
  virtual void
  operator()(const Eigen::Matrix<Scalar, n, 1>& x,
             Eigen::Matrix<Scalar, n, 1>* g) const = 0;
};

/// Dense Hessian functor

/// This functor evaluates the Hessian of the objective function
/// \f$ \nabla^2 f_0(\mathbf{x})\f$
/// Assumes the form:
/// \f[
/// \nabla^2 f_0(\mathbf{x}): \mathbb{R}^n \rightarrow \mathbb{R}^{n \times n}
/// \f]
template <typename Scalar, uint n>
struct HessianFunctor {
  /// Computes the Hessian
  virtual void
  operator()(const Eigen::Matrix<Scalar, n, 1>& x,
             Eigen::Matrix<Scalar, n, n>* h) const = 0;
};

/// Dummy constraint functor for unconstrained problems.

/// Helper object for uncosntrained problems
template <typename Scalar, uint n>
struct VoidConstraintFunctor : public ConstraintFunctor<Scalar, n> {
  /// Returns 0 always to comply with the ConstraintFunctor interface
  virtual Scalar operator()(const Eigen::Matrix<Scalar, n, 1>& x) const {
    return static_cast<Scalar>(0.0);
  }
};

/// Helper object for unconstrained problem

/// This helper object is a singleton object.
template <typename Scalar, uint m, uint n>
class NoConstraints : public Constraints<Scalar, m, n> {
 public:
  /// Returns the singleton of NoConstraints
  static NoConstraints<Scalar, m, n>& getInstance(void);

  // Non copyable
  // Copy constructor
  explicit NoConstraints(const NoConstraints<Scalar, m, n>& x) = delete;

  // Operator =
  NoConstraints<Scalar, m, n>&
  operator=(const NoConstraints<Scalar, m, n>& rhs) = delete;

  // Complying with interface
  virtual ConstraintFunctor<Scalar, n>& operator[](const unsigned int i) {
    return *c;
  }

 protected:
  static VoidConstraintFunctor<Scalar, n>* c;

 private:
  static void cleanUp(void);

  // Constructor
  NoConstraints(void) {
    c = new VoidConstraintFunctor<Scalar, n>;
    atexit(&cleanUp);
  }

  // Destructor
  ~NoConstraints(void) { }

  static NoConstraints<Scalar, m, n>* instance;
};

template <typename Scalar, uint m, uint n>
NoConstraints<Scalar, m, n>* NoConstraints<Scalar, m, n>::instance = nullptr;

template <typename Scalar, uint m, uint n>
VoidConstraintFunctor<Scalar, n>* NoConstraints<Scalar, m, n>::c = nullptr;

// Cleanup function
template <typename Scalar, uint m, uint n>
void NoConstraints<Scalar, m, n>::cleanUp(void) {
  delete c;
  c = nullptr;
}

template <typename Scalar, uint m, uint n>
NoConstraints<Scalar, m, n>& NoConstraints<Scalar, m, n>::getInstance(void) {
  if (!instance) {
    instance = new NoConstraints<Scalar, m, n>;
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
/// the problem fits in the stack. If it does not, then use the class
/// ProblemLS.
template <typename Scalar, uint m, uint n>
class Problem {
 public:
  const ObjectiveFunctor<Scalar, n>& objective; //< Objective functor
  const Constraints<Scalar, m, n>& constraints;  //< Constraints
  const GradientFunctor<Scalar, n>& gradient;  //< Gradient functor
  const HessianFunctor<Scalar, n>& hessian;  //< Hessian functor
  const Eigen::Matrix<Scalar, m, n> A;  //< Equality constraints matrix
  const Eigen::Matrix<Scalar, m, 1> b;  //< Equality constraints vector

  /// Constrained problem
  Problem(const ObjectiveFunctor<Scalar, n>& obj,
          const Constraints<Scalar, m, n>& con,
          const Eigen::Matrix<Scalar, m, n>& Aeq,
          const Eigen::Matrix<Scalar, m, 1>& beq,
          const GradientFunctor<Scalar, n>& g,
          const HessianFunctor<Scalar, n>& h) :
      objective(obj), constraints(con), gradient(g), hessian(h), A(Aeq), b(beq)
  { }

  /// Unconstrained problem
  Problem(const ObjectiveFunctor<Scalar, n>& obj,
          const GradientFunctor<Scalar, n>& g,
          const HessianFunctor<Scalar, n>& h) :
      objective(obj),
      constraints(NoConstraints<Scalar, m, n>::getInstance()),
      gradient(g), hessian(h) { }

  // Destructor
  virtual ~Problem(void) { }
};
}  // optimo
#endif  // OPTIMO_CORE_OBJECTS_H_
