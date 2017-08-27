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

#ifndef THEIA_SFM_TYPES_H_
#define THEIA_SFM_TYPES_H_

#include <Eigen/Core>
#include <cstdint>
#include <limits>
#include <stdint.h>
#include <utility>
#include <tuple>

namespace theia {

typedef uint32_t ViewId;
typedef uint32_t TrackId;
typedef uint32_t CameraIntrinsicsGroupId;
typedef std::pair<ViewId, ViewId> ViewIdPair;
typedef std::tuple<ViewId, ViewId, ViewId> ViewIdTriplet;

static const ViewId kInvalidViewId = std::numeric_limits<ViewId>::max();
static const TrackId kInvalidTrackId = std::numeric_limits<TrackId>::max();
static const CameraIntrinsicsGroupId kInvalidCameraIntrinsicsGroupId =
    std::numeric_limits<CameraIntrinsicsGroupId>::max();

// Used as the projection matrix type.
typedef Eigen::Matrix<double, 3, 4> Matrix3x4d;

}  // namespace theia

#endif  // THEIA_SFM_TYPES_H_
