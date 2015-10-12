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

#ifndef THEIA_IO_WRITE_MATCHES_DEPRECATED_H_
#define THEIA_IO_WRITE_MATCHES_DEPRECATED_H_

#include <string>
#include <vector>

#include "theia/matching/feature_matcher.h"
#include "theia/sfm/camera_intrinsics_prior.h"

namespace theia {

// Thi functions writes a Theia matches file. Matches files contain information
// needed to recreate a view graph and initialize a model before camera pose and
// structure estimation. This includes all two-view match information and
// optionally the two-view geometry as well. The format is as follows (in binary
// format with no newlines):
//
//   # views (uint32)
//   <View 1>
//   ...
//   <View n>
//   # image pair matches (uint64)
//   <Image pair match 1>
//   ...
//   <Image pair match k>
//
// Each view is comprised of the view name and metdata:
//   <Image name i> (uint32 for length of string, then the string)
//   image width, image height (2 ints)
//   <CameraIntrinsicsPrior focal length>
//   <CameraIntrinsicsPrior principal point[2]>
//
// The CameraIntrinsicsPrior written as <bool, double> corresponding to the
// is_set parameter and value parameter in //theia/sfm/camera_intrinsics_prior.h
//
// Each <Image pair match i> is formatted as follows:
//   <Image Index 1, Image Index 2> (2 uint32s)
//   <TwoViewInfo>
//   # features (uint32)
//   <Feature 1, Feature 2> (2 double[2])
//
// The <TwoViewInfo> is roughly the same format as the struct:
//   focal length 1 (double)
//   focal length 2 (double)
//   position_2 (double[3])
//   rotation_2 (double[3])
//   # 3D points (int)
//   # geometrically verified features (int)
//
// This file format is generic enough to allow for custom matching that can then
// be supplied to Theia for geometry estimation. If the custom matching does not
// produce two view geometry estimations then the TwoViewInfo can simply be
// ignored and reestimated using Theia.

// Writes the feature matches between view pairs as well as the two view
// geometry (i.e., TwoViewInfo) that describes the relative pose between the two
// views. The names of all views must be provided such that the image indices in
// the matches objects corresponds to the index of view_names. view_names should
// only store image names with extension and not the full image path.
// (e.g. abc.jpg and not /somepath/abc.jpg )
bool WriteMatchesAndGeometryDeprecated(
    const std::string& matches_file,
    const std::vector<std::string>& view_names,
    const std::vector<CameraIntrinsicsPrior>& camera_intrinsics_prior,
    const std::vector<ImagePairMatchDeprecated>& matches);

}  // namespace theia

#endif  // THEIA_IO_WRITE_MATCHES_DEPRECATED_H_
