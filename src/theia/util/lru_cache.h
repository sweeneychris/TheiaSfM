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

#ifndef THEIA_UTIL_LRU_CACHE_H_
#define THEIA_UTIL_LRU_CACHE_H_

#include <glog/logging.h>

#include <limits>
#include <list>
#include <mutex>
#include <unordered_map>
#include <utility>

#include "theia/util/map_util.h"
#include "theia/util/util.h"

namespace theia {

template <class KeyType, class ValueType>
class LRUCache {
  typedef std::list<KeyType> CacheList;
  typedef typename CacheList::iterator CacheListIterator;

 public:
  // Pass a function that performs the cache miss lookup (e.g., a read from
  // disk) that takes a key and returns a value.
  LRUCache(ValueType (*CacheMissLookup)(const KeyType&),
           const int max_cache_entries)
      : FetchEntryNotInCache(CacheMissLookup),
        max_cache_entries_(max_cache_entries) {
    CHECK_GT(max_cache_entries_, 0)
        << "The maximum number of cache entries must be greater than 0.";
    cache_misses_ = 0;
    cache_hits_ = 0;
  }

  // Fetch the entry and return the value. If the entry is in the cache then it
  // will be returned efficiently.
  virtual ValueType Fetch(const KeyType& key) {
    std::lock_guard<std::mutex> lock(mutex_);
    const auto it = cache_entries_map_.find(key);

    // If the value is not currently in the cache then we must fetch it.
    if (it == cache_entries_map_.end()) {
      ++cache_misses_;

      // Fetch the value for this key since it is not in the cache.
      const ValueType value = FetchEntryNotInCache(key);
      InsertIntoCache(key, value);
      return value;
    } else {
      ++cache_hits_;

      // If the entry was in the cache, we need to update the access record by
      // moving it to the back of the list.
      cache_entries_.splice(cache_entries_.end(),
                            cache_entries_,
                            it->second.second);

      return it->second.first;
    }
  }

  // Inserts a key-value pair into the cache, evicting the oldest entry if the
  // cache is at the maximum capacity. This method assumes that the key is not
  // already in the cache, and will CHECK-fail if the key already exists.
  virtual void Insert(const KeyType& key, const ValueType& value) {
    std::lock_guard<std::mutex> lock(mutex_);
    InsertIntoCache(key, value);
  }

  // Return if the key exists in the cache.
  virtual bool ExistsInCache(const KeyType& key) {
    return ContainsKey(cache_entries_map_, key);
  }

  // Varios statistics for the cache.
  int CacheCapacity() const { return max_cache_entries_; }
  int Size() const { return cache_entries_map_.size(); }
  int NumCacheMisses() const { return cache_misses_; }
  int NumCacheHits() const { return cache_hits_; }

 private:
  // Insert the key/value pair into the cache, evicting the oldest entry if
  // necessary.
  //
  // NOTE: This method is not thread-safe so any methods calling it must take
  // proper thread safety precautions.
  void InsertIntoCache(const KeyType& key, const ValueType& value) {
    // Ensure this method is only called on a cache miss.
    CHECK(!ContainsKey(cache_entries_map_, key));

    if (cache_entries_map_.size() == max_cache_entries_) {
      EvictOldestEntry();
    }

    // Insert the entry into the end of the accessor list (i.e. as the most
    // recently used).
    CacheListIterator it = cache_entries_.insert(cache_entries_.end(), key);

    // Add the entry to the map.
    cache_entries_map_.insert(std::make_pair(key, std::make_pair(value, it)));
  }

  // Evicts the oldest entry from the cache.
  //
  // NOTE: This method is not thread-safe so any methods calling it must take
  // proper thread safety precautions.
  void EvictOldestEntry() {
    // This method should never be called if the size of the cache is 0.
    CHECK_GT(cache_entries_.size(), 0);
    CHECK_GT(cache_entries_map_.size(), 0);

    const KeyType& evicted_key = *cache_entries_.begin();
    cache_entries_map_.erase(evicted_key);
    cache_entries_.pop_front();
  }

  // A function that takes in a KeyType as input and returns the ValueType. This
  // is utilized for cache misses and e.g., can implement a read from disk.
  ValueType (*FetchEntryNotInCache)(const KeyType&);

  // An ordered list that maintains the order in which cache entries have been
  // most recently added. The entries are oldest at the front and newest at the
  // back of the list.
  CacheList cache_entries_;

  // A lookup of keys to values. The lookup also provides an iterator directly
  // to the ordered cache list. This allows us to update the "least recently
  // used" quanitty in constant time.
  std::unordered_map<KeyType, std::pair<ValueType, CacheListIterator> >
      cache_entries_map_;

  // Maximum cache size.
  const int max_cache_entries_;

  // Some cache statistics.
  int cache_misses_, cache_hits_;

  std::mutex mutex_;

  DISALLOW_COPY_AND_ASSIGN(LRUCache);
};

}  // namespace theia

#endif  // THEIA_UTIL_LRU_CACHE_H_
