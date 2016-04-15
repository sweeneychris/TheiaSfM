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

#ifndef THEIA_UTIL_FILESYSTEM_H_
#define THEIA_UTIL_FILESYSTEM_H_

#include <string>
#include <vector>

namespace theia {

// Gets the filepath of all files matching the input wildcard. Returns true if
// the wildcard could be successfully evaluated and false otherwise (e.g. if the
// folder does not exist).
bool GetFilepathsFromWildcard(const std::string& filepath_with_wildcard,
                              std::vector<std::string>* filepaths);

// Extracts the filename from the filepath (i.e., removes all directory
// information). If with_extension is set to true then the extension is kept and
// output with the filename, otherwise the extension is removed.
bool GetFilenameFromFilepath(const std::string& filepath,
                             const bool with_extension,
                             std::string* filename);

// Returns the directory part of a given filepath.
bool GetDirectoryFromFilepath(const std::string& filepath,
                              std::string* directory);

// Returns true if the file exists, false otherwise.
bool FileExists(const std::string& filename);

// Returns true if the directory exists, false otherwise.
bool DirectoryExists(const std::string& directory);

// Creates the given directory.
bool CreateNewDirectory(const std::string& directory);

}  // namespace theia

#endif  // THEIA_UTIL_FILESYSTEM_H_
