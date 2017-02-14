// Copyright (C) 2017 The Regents of the University of California (Regents).
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
// Author: Aleksander Holynski (holynski@cs.washington.edu)

#include "theia/io/populate_image_sizes.h"

#include <glog/logging.h>

#include <cstdio>
#include <cstdlib>
#include <fstream>   // NOLINT
#include <iostream>  // NOLINT
#include <string>
#include <utility>
#include <vector>

#include "theia/image/image.h"
#include "theia/sfm/camera/camera.h"
#include "theia/sfm/reconstruction.h"
#include "theia/sfm/types.h"
#include "theia/sfm/view.h"
#include "theia/util/filesystem.h"
#include "theia/util/string.h"

namespace theia {

// Loads all images from the defined directory, and sets each of the
// recontruction's cameras to have an image size corresponding to the found
// image and a principal point at the center of that image.
bool PopulateImageSizesAndPrincipalPoints(const std::string& image_directory,
                                          Reconstruction* reconstruction) {
  CHECK_NOTNULL(reconstruction);
  std::string directory_with_slash = image_directory;
  AppendTrailingSlashIfNeeded(&directory_with_slash);
  const std::vector<ViewId> view_ids = reconstruction->ViewIds();
  for (int i = 0; i < view_ids.size(); i++) {
    const std::string file =
        directory_with_slash + reconstruction->View(view_ids[i])->Name();
    if (!FileExists(file)) {
      LOG(ERROR) << "Could not find " << file;
      return false;
    }
  }
  for (int i = 0; i < view_ids.size(); i++) {
    const std::string file =
        directory_with_slash + reconstruction->View(view_ids[i])->Name();
    const FloatImage image(file);
    CHECK_GT(image.Cols(), 0);
    Camera* camera = reconstruction->MutableView(view_ids[i])->MutableCamera();
    camera->SetImageSize(image.Cols(), image.Rows());
    camera->SetPrincipalPoint(image.Cols() / 2.0, image.Rows() / 2.0);
  }

  return true;
}

}  // namespace theia
