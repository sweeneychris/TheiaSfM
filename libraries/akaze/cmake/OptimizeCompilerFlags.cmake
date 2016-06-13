# Copyright (C) 2013 The Regents of the University of California (Regents).
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
#
#     * Neither the name of The Regents or University of California nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#
# Please contact the author of this library if you have any questions.
# Author: Chris Sweeney (cmsweeney@cs.ucsb.edu)
#
# These two methods are macros so that they can modify global parameters
# (functions cannot do that easily with CMake). Much of the compilation
# parameter setup was borrowed and possibly modified from Ceres Solver.
macro(OptimizeCompilerFlags)
  # Change the default build type from Debug to Release, while still
  # supporting overriding the build type.
  #
  # The CACHE STRING logic here and elsewhere is needed to force CMake
  # to pay attention to the value of these variables.
  if (NOT CMAKE_BUILD_TYPE)
    message("-- No build type specified; defaulting to CMAKE_BUILD_TYPE=Release.")
    set(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
  else (NOT CMAKE_BUILD_TYPE)
    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
      message("\n=================================================================================")
      message("\n-- Build type: Debug. Performance will be terrible!")
      message("-- Add -DCMAKE_BUILD_TYPE=Release to the CMake command line to get an optimized build.")
      message("\n=================================================================================")
    endif (CMAKE_BUILD_TYPE STREQUAL "Debug")
  endif (NOT CMAKE_BUILD_TYPE)

  # Set the default Akaze flags to an empty string.
  set (AKAZE_CXX_FLAGS)

  if (CMAKE_BUILD_TYPE STREQUAL "Release")
    if (CMAKE_COMPILER_IS_GNUCXX)
      # Linux
      if (CMAKE_SYSTEM_NAME MATCHES "Linux")
        if (NOT GCC_VERSION VERSION_LESS 4.2)
          set (AKAZE_CXX_FLAGS "${AKAZE_CXX_FLAGS} -march=native -mtune=native")
        endif (NOT GCC_VERSION VERSION_LESS 4.2)
      endif (CMAKE_SYSTEM_NAME MATCHES "Linux")
      # Mac OS X
      if (CMAKE_SYSTEM_NAME MATCHES "Darwin")
        set (AKAZE_CXX_FLAGS "${AKAZE_CXX_FLAGS} -msse3")
        # Use of -fast only applicable for Apple's GCC
        # Assume this is being used if GCC version < 4.3 on OSX
        execute_process(COMMAND ${CMAKE_C_COMPILER}
          ARGS ${CMAKE_CXX_COMPILER_ARG1} -dumpversion
          OUTPUT_VARIABLE GCC_VERSION
          OUTPUT_STRIP_TRAILING_WHITESPACE)
        if (GCC_VERSION VERSION_LESS 4.3)
          set (AKAZE_CXX_FLAGS "${AKAZE_CXX_FLAGS} -fast")
        endif (GCC_VERSION VERSION_LESS 4.3)
      endif (CMAKE_SYSTEM_NAME MATCHES "Darwin")
    endif (CMAKE_COMPILER_IS_GNUCXX)
    if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
      # Use of -O3 requires use of gold linker & LLVM-gold plugin, which might
      # well not be present / in use and without which files will compile, but
      # not link ('file not recognized') so explicitly check for support
      include(CheckCXXCompilerFlag)
      check_cxx_compiler_flag("-flto" HAVE_LTO_SUPPORT)
      if (HAVE_LTO_SUPPORT)
        message(STATUS "Enabling link-time optimization (-flto)")
        set(AKAZE_CXX_FLAGS "${AKAZE_CXX_FLAGS} -flto")
      else ()
        message(STATUS "Compiler/linker does not support link-time optimization (-flto), disabling.")
      endif (HAVE_LTO_SUPPORT)
    endif ()
  endif (CMAKE_BUILD_TYPE STREQUAL "Release")

  # Set c++ standard to c++11
  if (NOT MSVC)
    include(CheckCXXCompilerFlag)
    check_cxx_compiler_flag("-std=c++11" COMPILER_HAS_CXX11_FLAG)
    if (COMPILER_HAS_CXX11_FLAG)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
      if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
	# Mac OS X before Mavericks uses libstdc++ by default but does not support
	# c++11. Force it to use libc++.
	if (CMAKE_SYSTEM_NAME MATCHES "Darwin")
	  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
	endif (CMAKE_SYSTEM_NAME MATCHES "Darwin")
      endif (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    else (COMPILER_HAS_CXX11_FLAG)
      message(FATAL_ERROR "A compiler with C++11 support is required for Akaze.")
    endif (COMPILER_HAS_CXX11_FLAG)
  endif (NOT MSVC)

  set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${AKAZE_CXX_FLAGS}")

  # After the tweaks for the compile settings, disable some warnings on MSVC.
  if (MSVC)
    # Disable signed/unsigned int conversion warnings.
    add_definitions("/wd4018")
    # Disable warning about using struct/class for the same symobl.
    add_definitions("/wd4099")
    # Disable warning about the insecurity of using "std::copy".
    add_definitions("/wd4996")
    # Disable performance warning about int-to-bool conversion.
    add_definitions("/wd4800")
    # Disable performance warning about fopen insecurity.
    add_definitions("/wd4996")
    # Disable warning about int64 to int32 conversion. Disabling
    # this warning may not be correct; needs investigation.
    # TODO(keir): Investigate these warnings in more detail.
    add_definitions("/wd4244")
    # It's not possible to use STL types in DLL interfaces in a portable and
    # reliable way. However, that's what happens with Google Log and Google Flags
    # on Windows. MSVC gets upset about this and throws warnings that we can't do
    # much about. The real solution is to link static versions of Google Log and
    # Google Test, but that seems tricky on Windows. So, disable the warning.
    add_definitions("/wd4251")

    # Google Flags doesn't have their DLL import/export stuff set up correctly,
    # which results in linker warnings. This is irrelevant for Akaze, so ignore
    # the warnings.
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /ignore:4049")

    # Visual Studio has a limit to how many addresses and object can hold. This
    # can hobble templated classes that are large and result in compiler errors.
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj")

    # Update the C/CXX flags for MSVC to use either the static or shared
    # C-Run Time (CRT) library based on the user option: MSVC_USE_STATIC_CRT.
    list(APPEND C_CXX_FLAGS
      CMAKE_CXX_FLAGS
      CMAKE_CXX_FLAGS_DEBUG
      CMAKE_CXX_FLAGS_RELEASE
      CMAKE_CXX_FLAGS_MINSIZEREL
      CMAKE_CXX_FLAGS_RELWITHDEBINFO)

    foreach(FLAG_VAR ${C_CXX_FLAGS})
      if (MSVC_USE_STATIC_CRT)
        # Use static CRT.
        if (${FLAG_VAR} MATCHES "/MD")
          string(REGEX REPLACE "/MD" "/MT" ${FLAG_VAR} "${${FLAG_VAR}}")
        endif (${FLAG_VAR} MATCHES "/MD")
      else (MSVC_USE_STATIC_CRT)
        # Use shared, not static, CRT.
        if (${FLAG_VAR} MATCHES "/MT")
          string(REGEX REPLACE "/MT" "/MD" ${FLAG_VAR} "${${FLAG_VAR}}")
        endif (${FLAG_VAR} MATCHES "/MT")
      endif (MSVC_USE_STATIC_CRT)
    endforeach()

    # Tuple sizes of 10 are used by Gtest.
    add_definitions("-D_VARIADIC_MAX=10")
  endif (MSVC)

  if (UNIX)
    # GCC is not strict enough by default, so enable most of the warnings.
    set(CMAKE_CXX_FLAGS
      "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas -Wno-sign-compare -Wno-unused-parameter -Wno-missing-field-initializers")
  endif (UNIX)

endmacro(OptimizeCompilerFlags)
