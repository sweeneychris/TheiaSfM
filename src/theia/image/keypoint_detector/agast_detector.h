// Copyright (C) 2013 The Regents of the University of California (Regents).
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

#ifndef THEIA_IMAGE_KEYPOINT_DETECTOR_AGAST_DETECTOR_H_
#define THEIA_IMAGE_KEYPOINT_DETECTOR_AGAST_DETECTOR_H_

#include <agast/AstDetector.h>
#include <memory>
#include <vector>

#include "theia/image/keypoint_detector/keypoint.h"
#include "theia/image/keypoint_detector/keypoint_detector.h"
#include "theia/util/util.h"

namespace theia {
template<class T> class Image;
typedef Image<float> FloatImage;

// Detect keypoints using the AGAST method from "Adaptive and Generic Corner
// Detection Based on the Accelerated Segment Test" by Mair et. al.
class AgastDetector : public KeypointDetector {
 public:
  enum AstPattern {
    AGAST5_8,
    AGAST7_12D,
    AGAST7_12S,
    OAST9_16
  };

  // There are multiple patterns you can use for the AGAST detector. As the
  // pattern grows, the cost for computing keypoints increases. See the paper
  // for more details. The threshold should be specified, as well as whether
  // nonmaximum suppression should be run on the corners detected.
  AgastDetector(AstPattern pattern, int threshold, bool nonmax_suppression);
  AgastDetector(AstPattern pattern, int threshold)
      : AgastDetector(pattern, threshold, true) {}
  explicit AgastDetector(AstPattern pattern) : AgastDetector(pattern, 30) {}
  AgastDetector() : AgastDetector(AGAST5_8) {}

  ~AgastDetector() {}

  // Change the threshold for corner scores. You can do this without having to
  // create a new object, as this is performed per image.
  void SetThreshold(int threshold);

  // Detect the AGAST keypoints in the image.
  bool DetectKeypoints(const FloatImage& image,
                       std::vector<Keypoint>* keypoints);

 private:
  std::unique_ptr<agast::AstDetector> ast_detector_;

  bool nonmax_suppression_;
  DISALLOW_COPY_AND_ASSIGN(AgastDetector);
};

}  // namespace theia

#endif  // THEIA_IMAGE_KEYPOINT_DETECTOR_AGAST_DETECTOR_H_
