// The following code was taken from the OpenMVG library and is unmodified other
// than the namespace and header guard.
// To see the original source code, please visit OpenMVG's repository:
//     https://github.com/openMVG/openMVG
// The original software license is given below.

// ========================================================================== //
//
// Copyright (C) 2013 David Ok <david.ok8@gmail.com>
// Copyright (C) 2014 Pierre Moulon
//
// Adapted from DO++, a basic set of libraries in C++ for computer
// vision.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License v. 2.0. If a copy of the MPL was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// ========================================================================== //

#ifndef TIMER_HPP
#define TIMER_HPP

#ifdef HAVE_CXX11_CHRONO
#include <chrono>
#endif
#include <iostream>

namespace timer {

  //! \brief Timer class with microsecond accuracy.
  class Timer
  {
  public:
    //! Default constructor
    Timer();
    //! Reset the timer to zero.
    void reset();
    //! Returns the elapsed time in seconds.
    double elapsed() const;
    //! Returns the elapsed time in milliseconds.
    double elapsedMs() const;
  private:

#ifdef HAVE_CXX11_CHRONO
    std::chrono::high_resolution_clock::time_point start_;
#else
    double start_;
#ifdef _WIN32
    double frequency_;
#endif
#endif // HAVE_CXX11_CHRONO
  };

  // print the elapsed time
  std::ostream& operator << (std::ostream&, const Timer&);

} // namespace timer

#endif // TIMER_HPP
