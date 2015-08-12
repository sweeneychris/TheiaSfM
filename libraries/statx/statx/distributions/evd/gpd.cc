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

#include "statx/distributions/evd/gpd.h"
#include "statx/utils/ecdf.h"
#include "statx/distributions/evd/gpd_mle.h"
#include <vector>
#ifdef STATX_WITH_CERES
#include "statx/distributions/evd/gpd_ceres.h"
#endif

namespace statx {
namespace distributions {
namespace evd {

using std::vector;

// Main function calling either my MLE or CERES
// Parameters: sigma (scale) and xi (shape)
bool gpdfit(const vector<double>& data,
            double* xi,
            double* sigma,
            FitType fit_type) {
  if (data.empty()) return false;
  bool exit_flag = false;
  switch (fit_type) {
    default:
    case MLE:
      exit_flag = gpdfit_mle(data, xi, sigma);
      break;
    case QUANTILE_NLS:
#ifdef STATX_WITH_CERES
      exit_flag = gpdfit_ceres(data, xi, sigma);
#endif
      break;
  }
  return exit_flag;
}
}  // evd
}  // distributions
}  // statx
