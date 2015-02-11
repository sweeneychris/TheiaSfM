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
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVERCAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#include "statx/utils/ecdf.h"

#include <map>
#include <algorithm>

// TODO(vfragoso): Document me!
void statx::utils::ecdf(const std::vector<double>& y,
                        std::vector<double>* fx,
                        std::vector<double>* x) {
  if (y.empty() || !fx || !x) return;
  const size_t N = y.size();
  std::map<double, unsigned long> hist;
  size_t i;
  for (i = 0; i < N; i++) {
    if (hist.find(y[i]) == hist.end()) {
      hist[y[i]] = 0.0;
    }
    hist[y[i]]++;
  }
  fx->resize(hist.size());
  x->resize(hist.size());
  // S_hat(x) = Pi_{x_i < x} (1 - d_i/n_i)
  std::map<double, unsigned long>::const_iterator it;
  size_t d = 0;  // deaths (as in survival analysis)
  size_t r = N;  // samples at risk
  double prod = 1.0;
  i = 0;
  for (it = hist.begin(); it != hist.end(); ++it, i++) {
    d = it->second;  // deaths
    prod *= (1.0 - static_cast<double>(d)/static_cast<double>(r));  // s_hat
    r -= d;  // update number at risk
    (*x)[i] = it->first;  // x
    (*fx)[i] = 1.0 - prod;  // CDF = 1.0 - s_hat
  }
}
