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

#ifndef THEIA_SFM_TRANSFORMATION_ALIGN_RECONSTRUCTIONS_H_
#define THEIA_SFM_TRANSFORMATION_ALIGN_RECONSTRUCTIONS_H_

namespace theia {
class Reconstruction;

// Aligns the reconstructions so that their commons cameras have the closest
// positions in an L2 sense.
void AlignReconstructions(const Reconstruction& reconstruction1,
                          Reconstruction* reconstruction2);

// Aligns the reconstructions so that their commons cameras have the closest
// positions. This method is robust by using RANSAC to compute similarity
// transformations with inliers having a position distance less than
// robust_error_threshold. A final alignment is run on the inliers of the best
// RANSAC estimation.
void AlignReconstructionsRobust(
    const double robust_error_threshold,
    const Reconstruction& reconstruction1,
    Reconstruction* reconstruction2);

}  // namespace theia

#endif  // THEIA_SFM_TRANSFORMATION_ALIGN_RECONSTRUCTIONS_H_
