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

#ifndef THEIA_IMAGE_IMAGE_CACHE_H_
#define THEIA_IMAGE_IMAGE_CACHE_H_

#include <memory>
#include <string>

#include "theia/image/image.h"
#include "theia/util/lru_cache.h"

namespace theia {

// An LRU cache for retreiving images from disk. It is assumed that all images
// are held in the same directory and are referenced by their filename
// (including the extension).
class ImageCache {
public:
 // We require that all images are held in the same directory for
 // convenience. That directory must be specified here. and
 //
 // Images are held in an LRU cache so as to remain memory efficient for
 // out-of-core operations. This specifies the maximum number of images that
 // may be in the cache at any given time. This value may need to be adjusted
 // depending on how much memory is available.
 ImageCache(const std::string& image_directory,
            const int max_num_images_in_cache);
 ~ImageCache();

 // Returns the image corresponding to the view id, or a nullptr if the view
 // does not exist. It is up to the caller to determine if the image was
 // successfully fetched (i.e. if the returned image is not a nullptr). A
 // shared_ptr is used to monitor when images may be ejected from the cache and
 // prevent them from going out of scope while the caller is still using them.
 const std::shared_ptr<theia::FloatImage> FetchImage(
     const std::string& image_filename) const;

 const std::shared_ptr<theia::FloatImage> FetchGrayscaleImage(
     const std::string& image_filename) const;

 private:
  typedef LRUCache<std::string, std::shared_ptr<theia::FloatImage> >
      ImageLRUCache;

  // Method to fetch images from disk.
  std::shared_ptr<theia::FloatImage> FetchImagesFromDisk(
      const std::string& image_filename);

  // The directory where the images are stored.
  std::string image_directory_;

  // LRU Image cache.
  std::unique_ptr<ImageLRUCache> images_;
};

}  // namespace theia

#endif  // THEIA_IMAGE_IMAGE_CACHE_H_
