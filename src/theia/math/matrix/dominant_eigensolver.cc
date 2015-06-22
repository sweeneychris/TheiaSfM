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

#include "theia/math/matrix/dominant_eigensolver.h"

#include <Eigen/Core>
#include <glog/logging.h>
#include "theia/math/matrix/linear_operator.h"

namespace theia {

bool DominantEigensolver::Compute(double* eigenvalue,
                                  Eigen::VectorXd* eigenvector) const {
  Eigen::VectorXd previous_eigenvector(A_.Cols());
  previous_eigenvector.setRandom();
  for (int i = 0; i < options_.max_num_iterations; i++) {
    *eigenvector = previous_eigenvector.normalized();
    A_.RightMultiply(*eigenvector, &previous_eigenvector);
    *eigenvalue = eigenvector->dot(previous_eigenvector);

    const double error =
        (*eigenvector * (*eigenvalue) -  previous_eigenvector).stableNorm();
    VLOG(3) << "Iteration: " << i << "\tCurrent error = " << error;
    if (error < std::abs(*eigenvalue) * options_.tolerance) {
      VLOG(2) << "Power iterations converged after " << i + 1 << " iterations.";
      return true;
    }
  }

  return false;
}

}  // namespace theia
