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

#include "theia/image/descriptor/akaze_descriptor.h"

#include <Eigen/Core>
#include <algorithm>
#include <vector>

#include "akaze/src/AKAZE.h"
#include "glog/logging.h"

#include "theia/image/image.h"
#include "theia/image/descriptor/descriptor_extractor.h"
#include "theia/image/keypoint_detector/keypoint.h"

namespace theia {
bool AkazeDescriptorExtractor::ComputeDescriptor(const FloatImage& image,
                                                 const Keypoint& keypoint,
                                                 Eigen::VectorXf* descriptor) {
  LOG(FATAL) << "AKAZE must use its own Keypoints and so calling the AKAZE "
                "descriptor extractor with different keypoint cannot be used. "
                "Please use "
                "AkazeDescriptorExtractor::DetectAndExtractDescriptors() "
                "instead.";
}

bool AkazeDescriptorExtractor::DetectAndExtractDescriptors(
    const FloatImage& image,
    std::vector<Keypoint>* keypoints,
    std::vector<Eigen::VectorXf>* descriptors) {
  // Try to convert the image to grayscale and eigen type.
  const FloatImage& gray_image = image.AsGrayscaleImage();
  libAKAZE::RowMatrixXf img_32 =
      Eigen::Map<const Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic,
                                     Eigen::RowMajor> >(
          gray_image.Data(), gray_image.Rows(), gray_image.Cols());
  // REMOVE BELOW LINE (img_32 is already normalized)
  //img_32 /= 255.0;

  // Set the akaze options.
  libAKAZE::AKAZEOptions options;
  options.img_width = img_32.cols();
  options.img_height = img_32.rows();
  options.num_threads = 1;
  options.soffset = 1.6f;
  options.derivative_factor = 1.5f;
  options.omax = akaze_params_.maximum_octave_levels;
  options.nsublevels = akaze_params_.num_sublevels;
  options.dthreshold = akaze_params_.hessian_threshold;
  options.min_dthreshold = 0.00001f;

  options.diffusivity = libAKAZE::PM_G2;
  options.descriptor = libAKAZE::MSURF;
  options.descriptor_size = 0;
  options.descriptor_channels = 3;
  options.descriptor_pattern_size = 10;
  options.sderivatives = 1.0;

  options.kcontrast = 0.001f;
  options.kcontrast_percentile = 0.7f;
  options.kcontrast_nbins = 300;

  options.verbosity = false;

  // Extract features.
  std::vector<libAKAZE::AKAZEKeypoint> akaze_keypoints;
  libAKAZE::AKAZE evolution(options);
  evolution.Create_Nonlinear_Scale_Space(img_32);
  evolution.Feature_Detection(akaze_keypoints);

  // Compute descriptors.
  libAKAZE::AKAZEDescriptors akaze_descriptors;
  evolution.Compute_Descriptors(akaze_keypoints, akaze_descriptors);

  // Set the output keypoints.
  keypoints->reserve(akaze_keypoints.size());
  for (const auto& akaze_keypoint : akaze_keypoints) {
    Keypoint keypoint(akaze_keypoint.pt.x(), akaze_keypoint.pt.y(),
                      Keypoint::AKAZE);
    keypoint.set_scale(akaze_keypoint.size);
    keypoint.set_strength(akaze_keypoint.response);
    keypoint.set_orientation(akaze_keypoint.angle);
    keypoints->emplace_back(keypoint);
  }

  // Set the output descriptors.
  std::swap(*descriptors, akaze_descriptors.float_descriptor);
  return true;
}

}  // namespace theia
