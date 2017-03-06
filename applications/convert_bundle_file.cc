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

#include <Eigen/Core>
#include <glog/logging.h>
#include <gflags/gflags.h>
#include <theia/theia.h>

#include <algorithm>
#include <string>

DEFINE_string(lists_file, "", "Input bundle lists file.");
DEFINE_string(bundle_file, "", "Input bundle file.");
DEFINE_string(output_reconstruction_file, "",
              "Output reconstruction file in binary format.");
DEFINE_string(images_directory, "",
              "Directory of input images. This is used to extract the "
              "principal point and image dimensions since Bundler does not "
              "provide those.");
int main(int argc, char* argv[]) {
  google::InitGoogleLogging(argv[0]);
  THEIA_GFLAGS_NAMESPACE::ParseCommandLineFlags(&argc, &argv, true);

  // Load the reconstuction.
  theia::Reconstruction reconstruction;
  CHECK(theia::ReadBundlerFiles(FLAGS_lists_file,
                                FLAGS_bundle_file,
                                &reconstruction))
      << "Could not read Bundler files.";
  if (FLAGS_images_directory.size() > 0) {
    CHECK(theia::PopulateImageSizesAndPrincipalPoints(FLAGS_images_directory,
                                                      &reconstruction));
  } else {
    LOG(INFO) << "The image directory was not provided so the principal point "
                 "and image dimensions are assumed to be zero. Proceed with "
                 "caution!";
  }

  CHECK(WriteReconstruction(reconstruction, FLAGS_output_reconstruction_file))
      << "Could not write out reconstruction file.";
  return 0;
}
