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

#ifndef THEIA_UTIL_ENABLE_ENUM_BITMASK_OPERATORS_H_
#define THEIA_UTIL_ENABLE_ENUM_BITMASK_OPERATORS_H_

#include <type_traits>

namespace theia {

// In order to use bitmask operators with enum classes, you must implement the
// bitwise mask operator functions. We use enable_if to only enable this
// functionality if the field "enable" is set in the convenience struct to avoid
// enabling bitmask operators for every enum class.
//
// You may enable bitmask operators for an enum class with the
// ENABLE_ENUM_BITMASK_OPERATORS macro using code similar to the following:
//
//   enum class MyEnumClass {
//     val1 = 1,
//     val2 = 2,
//   };
//   ENABLE_ENUM_BITMASK_OPERATORS(MyEnumClass)

template <typename E>
struct EnableEnumBitmaskOperators {
    static const bool enable=false;
};

#define ENABLE_ENUM_BITMASK_OPERATORS(x) \
  template <>                            \
  struct EnableEnumBitmaskOperators<x> { \
    static const bool enable = true;     \
  };

template <typename E>
typename std::enable_if<EnableEnumBitmaskOperators<E>::enable, E>::type
operator|(E lhs, E rhs) {
  typedef typename std::underlying_type<E>::type underlying;
  return static_cast<E>(static_cast<underlying>(lhs) |
                        static_cast<underlying>(rhs));
}

template <typename E>
typename std::enable_if<EnableEnumBitmaskOperators<E>::enable, E>::type
operator&(E lhs, E rhs) {
  typedef typename std::underlying_type<E>::type underlying;
  return static_cast<E>(static_cast<underlying>(lhs) &
                        static_cast<underlying>(rhs));
}

template <typename E>
typename std::enable_if<EnableEnumBitmaskOperators<E>::enable, E>::type
operator^(E lhs, E rhs) {
  typedef typename std::underlying_type<E>::type underlying;
  return static_cast<E>(static_cast<underlying>(lhs) ^
                        static_cast<underlying>(rhs));
}

template <typename E>
typename std::enable_if<EnableEnumBitmaskOperators<E>::enable, E>::type
operator~(E lhs) {
  typedef typename std::underlying_type<E>::type underlying;
  return static_cast<E>(~static_cast<underlying>(lhs));
}

template <typename E>
typename std::enable_if<EnableEnumBitmaskOperators<E>::enable, E&>::type
operator|=(E& lhs, E rhs) {
  typedef typename std::underlying_type<E>::type underlying;
  lhs = static_cast<E>(static_cast<underlying>(lhs) |
                       static_cast<underlying>(rhs));
  return lhs;
}

template <typename E>
typename std::enable_if<EnableEnumBitmaskOperators<E>::enable, E&>::type
operator&=(E& lhs, E rhs) {
  typedef typename std::underlying_type<E>::type underlying;
  lhs = static_cast<E>(static_cast<underlying>(lhs) &
                       static_cast<underlying>(rhs));
  return lhs;
}

template <typename E>
typename std::enable_if<EnableEnumBitmaskOperators<E>::enable, E&>::type
operator^=(E& lhs, E rhs) {
  typedef typename std::underlying_type<E>::type underlying;
  lhs = static_cast<E>(static_cast<underlying>(lhs) ^
                       static_cast<underlying>(rhs));
  return lhs;
}

}  // namespace theia

#endif  // THEIA_UTIL_ENABLE_ENUM_BITMASK_OPERATORS_H_
