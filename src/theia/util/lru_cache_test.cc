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

#include "theia/util/lru_cache.h"

#include "gtest/gtest.h"

#include "theia/util/map_util.h"

namespace theia {

std::unordered_map<int, int> cache_lookup = {
    {0, 1}, {1, 47}, {2, 14}, {3, 101}, {4, 7}, {5, 29}};

// We will use this as our cache miss function for testing purposes.
int CacheMissLookup(const int& input) {
  return FindOrDie(cache_lookup, input);
}

TEST(LRUCache, Construtor) {
  const int kMaxCacheSize = 5;
  LRUCache<int, int> lru_cache(CacheMissLookup, kMaxCacheSize);
  EXPECT_EQ(lru_cache.CacheCapacity(), kMaxCacheSize);
  EXPECT_EQ(lru_cache.Size(), 0);
  EXPECT_EQ(lru_cache.NumCacheHits(), 0);
  EXPECT_EQ(lru_cache.NumCacheMisses(), 0);
}

TEST(LRUCache, FetchWorksWithOneEntry) {
  const int kMaxCacheSize = 5;
  LRUCache<int, int> lru_cache(CacheMissLookup, kMaxCacheSize);
  EXPECT_EQ(lru_cache.Fetch(0), FindOrDie(cache_lookup, 0));
  EXPECT_TRUE(lru_cache.ExistsInCache(0));
  // Ensure that our Fetch() works properly.
  EXPECT_EQ(lru_cache.Size(), 1);
  // Since the entries did not exist in the cache previously, it should be
  // registered as a cache miss when inserted.
  EXPECT_EQ(lru_cache.NumCacheMisses(), 1);
  EXPECT_EQ(lru_cache.NumCacheHits(), 0);

  // IF we fetch again, we should get the same results but with a cache hit.
  EXPECT_EQ(lru_cache.Fetch(0), FindOrDie(cache_lookup, 0));
  EXPECT_TRUE(lru_cache.ExistsInCache(0));
  // Ensure that our Fetch() works properly.
  EXPECT_EQ(lru_cache.Size(), 1);
  // Since the entries did not exist in the cache previously, it should be
  // registered as a cache miss when inserted.
  EXPECT_EQ(lru_cache.NumCacheMisses(), 1);
  EXPECT_EQ(lru_cache.NumCacheHits(), 1);
}

TEST(LRUCache, FetchWorksWithManyEntries) {
  const int kMaxCacheSize = 6;
  LRUCache<int, int> lru_cache(CacheMissLookup, kMaxCacheSize);

  // Add all elements to the cache.
  int num_elements_inserted = 0;
  for (const auto& entry : cache_lookup) {
    EXPECT_EQ(lru_cache.Fetch(entry.first), entry.second);
    EXPECT_TRUE(lru_cache.ExistsInCache(entry.first));
    ++num_elements_inserted;
    EXPECT_EQ(lru_cache.Size(), num_elements_inserted);
    // Since the entries did not exist in the cache previously, it should be
    // registered as a cache miss when inserted.
    EXPECT_EQ(lru_cache.NumCacheMisses(), num_elements_inserted);
    EXPECT_EQ(lru_cache.NumCacheHits(), 0);
  }

  // Attempt to fetch all elements and make sure they are cache hits.
  num_elements_inserted = 0;
  for (int i = 0; i < 5; i++) {
    for (const auto& entry : cache_lookup) {
      EXPECT_EQ(lru_cache.Fetch(entry.first), entry.second);
      EXPECT_TRUE(lru_cache.ExistsInCache(entry.first));
      ++num_elements_inserted;
      EXPECT_EQ(lru_cache.Size(), cache_lookup.size());
      // Since the entries did not exist in the cache previously, it should be
      // registered as a cache miss when inserted.
      EXPECT_EQ(lru_cache.NumCacheMisses(), cache_lookup.size());
      EXPECT_EQ(lru_cache.NumCacheHits(), num_elements_inserted);
    }
  }
}

TEST(LRUCache, FetchWithCacheMiss) {
  const int kMaxCacheSize = 1;
  LRUCache<int, int> lru_cache(CacheMissLookup, kMaxCacheSize);
  EXPECT_EQ(lru_cache.Fetch(0), FindOrDie(cache_lookup, 0));
  // IF we fetch again, we should evict the first entry.
  EXPECT_EQ(lru_cache.Fetch(1), FindOrDie(cache_lookup, 1));
  EXPECT_TRUE(lru_cache.ExistsInCache(1));
  EXPECT_FALSE(lru_cache.ExistsInCache(0));
  EXPECT_EQ(lru_cache.Size(), 1);
  EXPECT_EQ(lru_cache.NumCacheMisses(), 2);
  EXPECT_EQ(lru_cache.NumCacheHits(), 0);

  // Try to fetch element 0 again, it should be a miss.
  EXPECT_EQ(lru_cache.Fetch(0), FindOrDie(cache_lookup, 0));
  EXPECT_TRUE(lru_cache.ExistsInCache(0));
  EXPECT_FALSE(lru_cache.ExistsInCache(1));
  EXPECT_EQ(lru_cache.Size(), 1);
  EXPECT_EQ(lru_cache.NumCacheMisses(), 3);
  EXPECT_EQ(lru_cache.NumCacheHits(), 0);
}

}  // namespace theia
