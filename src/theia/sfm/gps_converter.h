// Copyright (C) 2016 The Regents of the University of California (Regents).
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

#ifndef THEIA_SFM_GPS_CONVERTER_H_
#define THEIA_SFM_GPS_CONVERTER_H_

#include <Eigen/Core>

namespace theia {

// This helper class contains static methods to convert between Geodetic
// coordinates (latitude, longitude, and altitude) and Earth-Center-Earth-Fixed
// (ECEF) coordinates. The method used is very fast, does not require special
// treatment at the poles or equator, and is extremely accurate. In tests by
// csweeney with randomly generated points, the maximum error of lat/lon is
// roughly 4e-16 and the maximum altitude error is roughly 3e-9 meters. These
// errors are smaller than the wavelength of visible light!!
//
// The original method was presenting in this paper:
// Olson, D.K. "Converting earth-Centered, Earth-Fixed Coordinates to Geodetic
// Coordinates," IEEE Transactions on Aerospace and Electronic Systems, Vol. 32,
// No. 1, January 1996, pp. 473-476.
class GPSConverter {
 public:
  GPSConverter() {}

  // Converts ECEF coordinates to GPS latitude, longitude, and altitude. ECEF
  // coordinates should be in meters. The returned latitude and longitude are in
  // degrees, and the altitude will be in meters.
  static Eigen::Vector3d ECEFToLLA(const Eigen::Vector3d& ecef);

  // Converts GPS latitude, longitude, and altitude to ECEF coordinates. The
  // latitude and longitude should be in degrees and the altitude in meters. The
  // returned ECEF coordinates will be in meters.
  static Eigen::Vector3d LLAToECEF(const Eigen::Vector3d& lla);
};

}  // namespace theia

#endif  // THEIA_SFM_GPS_CONVERTER_H_
