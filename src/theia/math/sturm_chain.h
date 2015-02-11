// Copyright (C) 2014 The Regents of the University of California (Regents).
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

#ifndef THEIA_MATH_STURM_CHAIN_H_
#define THEIA_MATH_STURM_CHAIN_H_

#include <Eigen/Core>
#include <vector>

namespace theia {

// A Sturm Chain (http://en.wikipedia.org/wiki/Sturm's_theorem) is a sequence of
// polynomials that reveals the GCD of polynomial p and its derivitive that is
// used to determine where the real polynomial roots lie.
class SturmChain {
 public:
  // Create a sturm chain from the polynomial.
  explicit SturmChain(const Eigen::VectorXd& polynomial);

  // Compute the number of sign changes when evaluating the sturm chain at
  // x. For an interval (a, b) the number of real roots in that interval is
  // NumSignChanges(a) - NumSignChanges(b).
  int NumSignChanges(const double x) const;

  // Compute the bounds of the real polynomial roots using Sameulson's
  // inequality.
  void ComputeRootBounds(double* lower_bound, double* upper_bound);

  // The difference in the number of sign changes at infinity tell us the
  // theoretical number of real roots.
  int NumSignChangesAtInfinity() const;
  int NumSignChangesAtNegativeInfinity() const;

 private:
  std::vector<Eigen::VectorXd> sturm_chain_;
};

}  // namespace theia

#endif  // THEIA_MATH_STURM_CHAIN_H_
