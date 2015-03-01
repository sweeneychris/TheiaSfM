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

#ifndef THEIA_IO_READ_1DSFM_H_
#define THEIA_IO_READ_1DSFM_H_

#include <string>

namespace theia {

class Reconstruction;
class ViewGraph;

// Reads a 1dSfM dataset from http://www.cs.cornell.edu/projects/1dsfm/ and the
// publication: "Robust Global Translations with 1DSfM" by Kyle Wilson and Noah
// Snavely (ECCV 2014). This dataset consists of initial two view matches and
// geometry between images as well as tracks. This geometric information can be
// used to initialize a 3D reconstruction and, as such, can serve as a benchmark
// tool for structure from motion.
//
// The dataset directory should *not* contain a trailing slash and must contain
// 5 files. These files are included with the 1dSfM datasets by default:
//   lists.txt: The Bundler lists files containing image names and possible EXIF
//     calibration.
//   EGs.txt: The epipolar geometries between matching images.
//   cc.txt: The largest connected component of the model.
//   coords.txt: The image feature coordinates.
//   tracks.txt: The tracks that are created from verified image matches.
//
// Given these inputs, a Reconstruction and ViewGraph will be initialized to
// exactly correspond to the input data. These objects may then be passed to the
// ReconstructionEstimator to estimate global poses and triangulate 3D points.
//
// Returns true on success, and false if one or more of the input files could
// not be found.
bool Read1DSFM(const std::string& dataset_directory,
               Reconstruction* reconstruction,
               ViewGraph* view_graph);

}  // namespace theia

#endif  // THEIA_IO_READ_1DSFM_H_
