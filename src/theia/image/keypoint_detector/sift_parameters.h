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
// Author: Victor Fragoso (vfragoso@cs.ucsb.edu)

#ifndef THEIA_IMAGE_KEYPOINT_DETECTOR_SIFT_PARAMETERS_H_
#define THEIA_IMAGE_KEYPOINT_DETECTOR_SIFT_PARAMETERS_H_

// Sift blob feature detector parameters. Since the Sift implementation is based
// on the VLFeat one, please visit (http://www.vlfeat.org/api/sift.html) for
// getting more info about the parameters.
struct SiftParameters {
  SiftParameters(int num_octaves,
                 int num_levels,
                 int first_octave,
                 float edge_threshold,
                 float peak_threshold) :
      num_octaves(num_octaves), num_levels(num_levels),
      first_octave(first_octave), edge_threshold(edge_threshold),
      peak_threshold(peak_threshold) {}

  SiftParameters(int num_octaves, int num_levels, int first_octave) :
      SiftParameters(num_octaves, num_levels, first_octave,
                     5.0f, 0.5f) {}
  ~SiftParameters() {}

  // Parameters.
  // Blob feature detector params.
  int num_octaves = -1;
  int num_levels = 3;
  int first_octave = 0;
  float edge_threshold = 5.0f;
  float peak_threshold = 0.5f;
  // Descriptor parameters.
  bool root_sift = true;
};

#endif  // THEIA_IMAGE_KEYPOINT_DETECTOR_SIFT_PARAMETERS_H_
