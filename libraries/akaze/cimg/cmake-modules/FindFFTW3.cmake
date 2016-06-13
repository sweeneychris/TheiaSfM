# - Try to find FFTW3.
# Usage: find_package(FFTW3 [COMPONENTS [single double long-double threads]])
#
# Variables used by this module:
#  FFTW3_ROOT_DIR             - FFTW3 root directory
# Variables defined by this module:
#  FFTW3_FOUND                - system has FFTW3
#  FFTW3_INCLUDE_DIR          - the FFTW3 include directory (cached)
#  FFTW3_INCLUDE_DIRS         - the FFTW3 include directories
#                               (identical to FFTW3_INCLUDE_DIR)
#  FFTW3[FL]?_LIBRARY         - the FFTW3 library - double, single(F), 
#                               long-double(L) precision (cached)
#  FFTW3[FL]?_THREADS_LIBRARY - the threaded FFTW3 library - double, single(F), 
#                               long-double(L) precision (cached)
#  FFTW3_LIBRARIES            - list of all FFTW3 libraries found

# Copyright (C) 2009-2010
# ASTRON (Netherlands Institute for Radio Astronomy)
# P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
#
# This file is part of the LOFAR software suite.
# The LOFAR software suite is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# The LOFAR software suite is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with the LOFAR software suite. If not, see <http://www.gnu.org/licenses/>.
#
# $Id: FindFFTW3.cmake 15918 2010-06-25 11:12:42Z loose $

# Use double precision by default.
if(FFTW3_FIND_COMPONENTS MATCHES "^$")
  set(_components double)
else()
  set(_components ${FFTW3_FIND_COMPONENTS})
endif()

# Loop over each component.
set(_libraries)
foreach(_comp ${_components})
  if(_comp STREQUAL "single")
    list(APPEND _libraries fftw3f)
  elseif(_comp STREQUAL "double")
    list(APPEND _libraries fftw3)
  elseif(_comp STREQUAL "long-double")
    list(APPEND _libraries fftw3l)
  elseif(_comp STREQUAL "threads")
    set(_use_threads ON)
  else(_comp STREQUAL "single")
    message(FATAL_ERROR "FindFFTW3: unknown component `${_comp}' specified. "
      "Valid components are `single', `double', `long-double', and `threads'.")
  endif(_comp STREQUAL "single")
endforeach(_comp ${_components})

# If using threads, we need to link against threaded libraries as well.
if(_use_threads)
  set(_thread_libs)
  foreach(_lib ${_libraries})
    list(APPEND _thread_libs ${_lib}_threads)
  endforeach(_lib ${_libraries})
  set(_libraries ${_thread_libs} ${_libraries})
endif(_use_threads)

# Keep a list of variable names that we need to pass on to
# find_package_handle_standard_args().
set(_check_list)

# Search for all requested libraries.
foreach(_lib ${_libraries})
  string(TOUPPER ${_lib} _LIB)
  find_library(${_LIB}_LIBRARY ${_lib}
    HINTS ${FFTW3_ROOT_DIR} PATH_SUFFIXES lib)
  mark_as_advanced(${_LIB}_LIBRARY)
  list(APPEND FFTW3_LIBRARIES ${${_LIB}_LIBRARY})
  list(APPEND _check_list ${_LIB}_LIBRARY)
endforeach(_lib ${_libraries})

# Search for the header file.
