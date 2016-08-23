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

#include "theia/image/image_cache.h"

#include <memory>
#include <string>

#include "theia/image/image.h"
#include "theia/util/lru_cache.h"
#include "theia/util/string.h"

namespace theia {

ImageCache::ImageCache(const std::string& image_directory,
                       const int max_num_images_in_cache)
    : image_directory_(image_directory) {
  AppendTrailingSlashIfNeeded(&image_directory_);

  // Set up the LRU image cache.
  const std::function<std::shared_ptr<FloatImage>(const std::string)>
      fetch_images = std::bind(&ImageCache::FetchImagesFromDisk,
                               this,
                               std::placeholders::_1);
  images_.reset(new ImageLRUCache(fetch_images, max_num_images_in_cache));
}

ImageCache::~ImageCache() {}

const std::shared_ptr<FloatImage> ImageCache::FetchImage(
    const std::string& image_filename) const {
  const std::shared_ptr<FloatImage> image = images_->Fetch(image_filename);
  return image;
}

std::shared_ptr<FloatImage> ImageCache::FetchImagesFromDisk(
    const std::string& image_filename) {
  const std::string image_filepath = image_directory_ + image_filename;
  std::shared_ptr<theia::FloatImage> image =
      std::make_shared<FloatImage>(image_filepath);
  return image;
}

}  // namespace theia
