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

#ifndef OPTIMO_SOLVERS_SOLVER_H_
#define OPTIMO_SOLVERS_SOLVER_H_

namespace optimo {
namespace solvers {

typedef unsigned int uint;

/// This enum captures the possible termination states that a solver can have.
enum TERMINATION_TYPE {
  SOLVED = 0,  //< Solver succesfully solved the problem
  FAIL_NAN_INF = 1,  //< Solver encountered NaN/Inf
  MAX_ITERATIONS = 2,  //< Solver reached max. number of iterations
  NOT_SOLVED = 3,  //< Solver could not solve problem
  // Params used to map the SuiteSparse return types for sparse-matrix
  // based solvers.
  INVALID_ARGUMENTS = 4,  //< Solver received invalid arguments/parameters
  NUMERICAL_ISSUE = 5,  //< Solver encountered a numerical issue
  NO_CONVERGENCE = 6,  //< Solver did not converge
  INFEASIBLE_STARTING_POINT = 7  //< Starting point is infeasible
};

/// Abstract class for a Solver algorithm.
template <typename Scalar>
class Solver {
 public:
  /// Various parameters for line search and tolerances.
  struct Options {
    // Line search params (backtracking params)
    Scalar alpha_ = 0.01;  ///< Backtracking (Line search) parameter
    Scalar beta_ = 0.5;  ///< Backtracking (line search) parameter
    // Primal dual and Newton params
    Scalar mu_ = 10;  ///< For primal dual computation
    Scalar eps_feas_ = 1e-8;  ///< Tolerance for primal dual residuals
    Scalar epsilon_ = 1e-6;  ///< Tolerance
    uint max_iter_ = 1000;  ///< Maximum number of iterations
  } options; ///< Solver parameters
};
}  // solvers
}  // optimo
#endif  // OPTIMO_SOLVERS_SOLVER_H_
