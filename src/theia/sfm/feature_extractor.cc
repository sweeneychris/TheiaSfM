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

#include "theia/sfm/feature_extractor.h"

#include <Eigen/Core>
#include <string>
#include <vector>

#include "theia/image/descriptor/brief_descriptor.h"
#include "theia/image/descriptor/brisk_descriptor.h"
#include "theia/image/descriptor/descriptor_extractor.h"
#include "theia/image/descriptor/freak_descriptor.h"
#include "theia/image/descriptor/sift_descriptor.h"
#include "theia/image/image.h"

namespace theia {

std::unique_ptr<DescriptorExtractor>
FeatureExtractor::CreateDescriptorExtractor(
    const DescriptorExtractorType& descriptor_extractor_type,
    const SiftParameters& sift_parameters) {
  std::unique_ptr<DescriptorExtractor> descriptor_extractor;
  switch (descriptor_extractor_type) {
    case DescriptorExtractorType::SIFT:
      descriptor_extractor.reset(new SiftDescriptorExtractor(sift_parameters));
      break;
    case DescriptorExtractorType::BRIEF:
      descriptor_extractor.reset(new BriefDescriptorExtractor);
      break;
    case DescriptorExtractorType::BRISK:
      descriptor_extractor.reset(new BriskDescriptorExtractor);
      break;
    case DescriptorExtractorType::FREAK:
      descriptor_extractor.reset(new FreakDescriptorExtractor);
      break;
    default:
      LOG(ERROR) << "Invalid Descriptor Extractor specified.";
  }
  CHECK(descriptor_extractor->Initialize())
      << "Could not initialize the Descriptor Extractor";
  return descriptor_extractor;
}

}  // namespace theia
