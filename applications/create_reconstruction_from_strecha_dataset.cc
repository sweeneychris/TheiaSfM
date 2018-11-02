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
//
// This program is a useful tool for converting camera files from the Strecha
// Multi-View Stereo dataset into Theia reconstructions. The output
// reconstruction may serve as ground truth for experiments that utilize the
// Theia library. The datasets may be found at the link below, and the relevant
// files for this program are the "camera" files.
//
// http://cvlabwww.epfl.ch/data/multiview/denseMVS.html

#include <Eigen/Core>
#include <glog/logging.h>
#include <gflags/gflags.h>
#include <theia/theia.h>

#include <fstream>  // NOLINT
#include <string>

DEFINE_string(input_dataset_directory, "",
              "Directory containing camera files from the Strecha MVS dataset. "
              "Do not include a trailing slash.");
DEFINE_string(output_reconstruction, "",
              "Filepath of the output Theia reconstruction file generated.");

int main(int argc, char* argv[]) {
  google::InitGoogleLogging(argv[0]);
  THEIA_GFLAGS_NAMESPACE::ParseCommandLineFlags(&argc, &argv, true);

  // Add each camera in the Strecha MVS dataset to the reconstruciton.
  theia::Reconstruction reconstruction;
  CHECK(theia::ReadStrechaDataset(FLAGS_input_dataset_directory, &reconstruction));
  CHECK(theia::WriteReconstruction(reconstruction, FLAGS_output_reconstruction))
      << "Could not write the reconstruction file to: "
      << FLAGS_output_reconstruction;

  return 0;
}
