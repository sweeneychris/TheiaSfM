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

#include <glog/logging.h>
#include <gflags/gflags.h>
#include <theia/theia.h>
#include <fstream>  // NOLINT
#include <string>
#include <vector>

// Input/output files.
DEFINE_string(images, "", "Wildcard of images to reconstruct.");
DEFINE_string(output_calibration_file, "",
              "Calibration file containing image calibration data.");
DEFINE_bool(initialize_uncalibrated_images_with_median_viewing_angle, true,
            "Images with no EXIF information initialize the focal length based "
            "on a focal length corresponding to a median viewing angle.");

int main(int argc, char *argv[]) {
  THEIA_GFLAGS_NAMESPACE::ParseCommandLineFlags(&argc, &argv, true);
  google::InitGoogleLogging(argv[0]);

  std::vector<std::string> image_files;
  CHECK(theia::GetFilepathsFromWildcard(FLAGS_images, &image_files))
      << "Could not find images that matched the filepath: " << FLAGS_images
      << ". NOTE that the ~ filepath is not supported.";

  std::ofstream ofs(FLAGS_output_calibration_file, std::ios::out);
  if (!ofs.is_open()) {
    LOG(ERROR) << "Could not open the calibration file: "
               << FLAGS_output_calibration_file << " for writing.";
    return false;
  }

  theia::ExifReader exif_reader;

  //   image_name focal_length ppx ppy aspect_ratio skew k1 k2
  for (int i = 0; i < image_files.size(); i++) {
    std::string image_name;
    theia::GetFilenameFromFilepath(image_files[i], true, &image_name);

    theia::CameraIntrinsicsPrior prior;
    CHECK(exif_reader.ExtractEXIFMetadata(image_files[i], &prior))
        << "Could not open " << image_files[i] << " for reading.";

    // Only write the calibration for images with a focal length that was
    // extracted.
    if (!prior.focal_length.is_set) {
      if (FLAGS_initialize_uncalibrated_images_with_median_viewing_angle) {
        // Set the focal length based on a median viewing angle.
        prior.focal_length.is_set = true;
        prior.focal_length.value[0] =
            1.2 * static_cast<double>(
                      std::max(prior.image_width, prior.image_height));
      } else {
        LOG(INFO) << image_name << " did not contain an EXIF focal length.";
        continue;
      }
    }

    // We write the default values for aspect ratio, skew, and radial distortion
    // since those cannot be recovered from EXIF.
    LOG(INFO) << image_name << " has an EXIF focal length of "
              << prior.focal_length.value[0];
    ofs << image_name << " " << prior.focal_length.value[0] << " "
        << prior.principal_point.value[0] << " "
        << prior.principal_point.value[1] << " 1.0 0.0 0.0 0.0\n";
  }
  ofs.close();
}
