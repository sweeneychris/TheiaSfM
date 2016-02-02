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

DEFINE_string(reconstruction, "", "Theia Reconstruction file.");
DEFINE_string(images, "",
              "A filepath wildcard specifying all images that were used in the "
              "reconstruction.");
DEFINE_string(pmvs_working_directory, "",
              "A directory to store the necessary pmvs files.");
DEFINE_int32(num_threads, 1, "Number of threads used in PMVS.");

void CreateDirectoryIfDoesNotExist(const std::string& directory) {
  if (!theia::DirectoryExists(directory)) {
    CHECK(theia::CreateDirectory(directory))
        << "Could not create the directory: " << directory;
  }
}

int WriteCamerasToPMVS(const theia::Reconstruction& reconstruction) {
  const std::string txt_dir = FLAGS_pmvs_working_directory + "/txt";
  CreateDirectoryIfDoesNotExist(txt_dir);
  const std::string visualize_dir = FLAGS_pmvs_working_directory + "/visualize";

  std::vector<std::string> image_files;
  CHECK(theia::GetFilepathsFromWildcard(FLAGS_images, &image_files))
      << "Could not find images that matched the filepath: " << FLAGS_images
      << ". NOTE that the ~ filepath is not supported.";
  CHECK_GT(image_files.size(), 0) << "No images found in: " << FLAGS_images;

  // Format for printing eigen matrices.
  const Eigen::IOFormat unaligned(Eigen::StreamPrecision, Eigen::DontAlignCols);

  int current_image_index = 0;
  for (int i = 0; i < image_files.size(); i++) {
    std::string image_name;
    CHECK(theia::GetFilenameFromFilepath(image_files[i], true, &image_name));
    const theia::ViewId view_id = reconstruction.ViewIdFromName(image_name);
    if (view_id == theia::kInvalidViewId) {
      continue;
    }

    LOG(INFO) << "Exporting parameters for image: " << image_name;

    // Copy the image into a jpeg format with the filename in the form of
    // %08d.jpg.
    const std::string new_image_file = theia::StringPrintf(
        "%s/%08d.jpg", visualize_dir.c_str(), current_image_index);
    theia::Image<float> old_image(image_files[i]);
    old_image.Write(new_image_file);

    // Write the camera projection matrix.
    const std::string txt_file = theia::StringPrintf(
        "%s/%08d.txt", txt_dir.c_str(), current_image_index);
    const theia::Camera camera = reconstruction.View(view_id)->Camera();

    theia::Matrix3x4d projection_matrix;
    camera.GetProjectionMatrix(&projection_matrix);
    std::ofstream ofs(txt_file);
    ofs << "CONTOUR" << std::endl;
    ofs << projection_matrix.format(unaligned);
    ofs.close();

    ++current_image_index;
  }

  return current_image_index;
}

void WritePMVSOptions(const std::string& working_dir,
                      const int num_images) {
  std::ofstream ofs(working_dir + "/pmvs_options.txt");
  ofs << "level 1" << std::endl;
  ofs << "csize 2" << std::endl;
  ofs << "threshold 0.7" << std::endl;
  ofs << "wsize 7" << std::endl;
  ofs << "minImageNum 3" << std::endl;
  ofs << "CPU " << FLAGS_num_threads << std::endl;
  ofs << "setEdge 0" << std::endl;
  ofs << "useBound 0" << std::endl;
  ofs << "useVisData 0" << std::endl;
  ofs << "sequence -1" << std::endl;
  ofs << "timages -1 0 " << num_images << std::endl;
  ofs << "oimages 0" << std::endl;
}

int main(int argc, char* argv[]) {
  google::InitGoogleLogging(argv[0]);
  THEIA_GFLAGS_NAMESPACE::ParseCommandLineFlags(&argc, &argv, true);

  // Set up output directories.
  CreateDirectoryIfDoesNotExist(FLAGS_pmvs_working_directory);
  const std::string visualize_dir = FLAGS_pmvs_working_directory + "/visualize";
  CreateDirectoryIfDoesNotExist(visualize_dir);
  const std::string txt_dir = FLAGS_pmvs_working_directory + "/txt";
  CreateDirectoryIfDoesNotExist(txt_dir);
  const std::string models_dir = FLAGS_pmvs_working_directory + "/models";
  CreateDirectoryIfDoesNotExist(models_dir);

  theia::Reconstruction reconstruction;
  CHECK(theia::ReadReconstruction(FLAGS_reconstruction, &reconstruction))
      << "Could not read Reconstruction files.";

  const int num_cameras = WriteCamerasToPMVS(reconstruction);
  WritePMVSOptions(FLAGS_pmvs_working_directory, num_cameras);

  const std::string lists_file = FLAGS_pmvs_working_directory + "/list.txt";
  const std::string bundle_file =
      FLAGS_pmvs_working_directory + "/bundle.rd.out";
  CHECK(theia::WriteBundlerFiles(reconstruction, lists_file, bundle_file));

  return 0;
}
