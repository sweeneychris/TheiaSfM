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

#include <Eigen/Core>
#include <glog/logging.h>
#include <gflags/gflags.h>
#include <theia/theia.h>

#include <algorithm>
#include <memory>
#include <string>

DEFINE_string(input_reconstruction, "",
              "Input reconstruction file with distorted cameras");
DEFINE_string(output_reconstruction, "",
              "Output reconstruction file with undistorted cameras");
DEFINE_string(input_image_directory, "",
              "Directory containing the input distorted images.");
DEFINE_string(output_image_directory, "",
              "Directory containing the output undistorted images.");
DEFINE_int32(num_threads, 1, "Number of threads to use for undistortion.");

void UndistortImageAndWriteToFile(const std::string input_image_directory,
                                  const std::string output_image_directory,
                                  const theia::ViewId view_id,
                                  theia::Reconstruction* reconstruction) {
  theia::View* view = reconstruction->MutableView(view_id);
  LOG(INFO) << "Undistorting image " << view->Name();

  // Skip this image if the view is unestimated or does not exist.
  if (view == nullptr || !view->IsEstimated()) {
    return;
  }

  // Undistort the image.
  const theia::FloatImage distorted_image(input_image_directory + view->Name());
  theia::View undistorted_view;
  theia::FloatImage undistorted_image;
  CHECK(theia::UndistortView(*view,
                             distorted_image,
                             &undistorted_view,
                             &undistorted_image))
      << "Could not undistort image: " << view->Name();

  // Update the view to the undistorted view.
  *view = undistorted_view;

  // Save the image to the output directory.
  LOG(INFO) << "Writing image " << view->Name();
  undistorted_image.Write(output_image_directory + view->Name());
}

int main(int argc, char* argv[]) {
  google::InitGoogleLogging(argv[0]);
  THEIA_GFLAGS_NAMESPACE::ParseCommandLineFlags(&argc, &argv, true);

  // Load the reconstruction.
  theia::Reconstruction reconstruction;
  CHECK(theia::ReadReconstruction(FLAGS_input_reconstruction, &reconstruction))
      << "Could not read reconstruction file.";

  // Get the image input and output directories.
  std::string input_image_directory = FLAGS_input_image_directory;
  theia::AppendTrailingSlashIfNeeded(&input_image_directory);
  std::string output_image_directory = FLAGS_output_image_directory;
  theia::AppendTrailingSlashIfNeeded(&output_image_directory);

  // Undistort images in parallel.
  theia::ThreadPool pool(FLAGS_num_threads);
  const auto& view_ids = reconstruction.ViewIds();
  for (const theia::ViewId view_id : view_ids) {
    pool.Add(UndistortImageAndWriteToFile,
             input_image_directory,
             output_image_directory,
             view_id,
             &reconstruction);
  }
  pool.WaitForTasksToFinish();

  theia::WriteReconstruction(reconstruction, FLAGS_output_reconstruction);
  return 0;
}
