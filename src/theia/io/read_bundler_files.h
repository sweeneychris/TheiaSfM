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

#ifndef THEIA_IO_READ_BUNDLER_FILES_H_
#define THEIA_IO_READ_BUNDLER_FILES_H_

#include <string>

namespace theia {

class Reconstruction;

// Loads all information from a bundler file and the corresponding lists file
// that provide the image names and EXIF focal lengths. The bundler file
// includes 3D points, camera poses, camera intrinsics, descriptors, and 2D-3D
// matches. This method will not load the descriptors from the sift key files
// and it will ignore any SIFT descriptor information.
//
// Input params are as follows:
//   lists_file: the list of image names from the bundler list file. The focal
//       length is provided from EXIF data or set to 0 if no EXIF focal length
//       is available.
//   bundler_file: the file output by bundler containing the 3D reconstruction,
//       camera poses, feature locations, and correspondences. Usually the file
//       is named bundle.out
//   reconstruction: A Theia Reconstruction containing the camera, track, and
//       point cloud information. See theia/sfm/reconstruction.h for more
//       information.
bool ReadBundlerFiles(const std::string& lists_file,
                      const std::string& bundle_file,
                      Reconstruction* reconstruction);

// Bundler files & image lists don't usually contain image sizes.
// This function loads the images with names defined in the reconstruction from
// the 'image_directory' folder.
// If any of the files defined in the reconstruction do not exist, the function
// will return false (and no values will be changed in the reconstruction),
// otherwise the function will return true.
// This function is to be called after ReadBundlerFiles().
// Assumes principal points to be at the image center.
//
// Input params are as follows:
//   image_directory: The directory containing all the image files from the
//   reconstruction.
//   reconstruction: A Theia Reconstruction containing the camera, track, and
//       point cloud information. See theia/sfm/reconstruction.h for more
//       information.
bool PopulateImageSizesAndPrincipalPoints(const std::string& image_directory,
                                          Reconstruction* reconstruction);

}  // namespace theia

#endif  // THEIA_IO_READ_BUNDLER_FILES_H_
