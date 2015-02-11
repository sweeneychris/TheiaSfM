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
# Much of this file was modified from Ceres Solver which has the license below:
#
# Ceres Solver - A fast non-linear least squares minimizer
# Copyright 2010, 2011, 2012 Google Inc. All rights reserved.
# http://code.google.com/p/ceres-solver/
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright notice,
#   this list of conditions and the following disclaimer.
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
# * Neither the name of Google Inc. nor the names of its contributors may be
#   used to endorse or promote products derived from this software without
#   specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#
# Author: keir@google.com (Keir Mierle)

# These two methods are macros so that they can modify global parameters
# (functions cannot do that easily with CMake). Much of the compilation
# parameter setup was borrowed and possibly modified from Ceres Solver.
macro(OptimizeTheiaCompilerFlags)
  # Change the default build type from Debug to Release, while still
  # supporting overriding the build type.
  #
  # The CACHE STRING logic here and elsewhere is needed to force CMake
  # to pay attention to the value of these variables.
  IF (NOT CMAKE_BUILD_TYPE)
    MESSAGE("-- No build type specified; defaulting to CMAKE_BUILD_TYPE=Release.")
    SET(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
  ELSE (NOT CMAKE_BUILD_TYPE)
    IF (CMAKE_BUILD_TYPE STREQUAL "Debug")
      MESSAGE("\n=================================================================================")
      MESSAGE("\n-- Build type: Debug. Performance will be terrible!")
      MESSAGE("-- Add -DCMAKE_BUILD_TYPE=Release to the CMake command line to get an optimized build.")
      MESSAGE("\n=================================================================================")
    ENDIF (CMAKE_BUILD_TYPE STREQUAL "Debug")
  ENDIF (NOT CMAKE_BUILD_TYPE)

  # Set the default Theia flags to an empty string.
  SET (THEIA_CXX_FLAGS)

  IF (CMAKE_BUILD_TYPE STREQUAL "Release")
    IF (CMAKE_COMPILER_IS_GNUCXX)
      # Linux
      IF (CMAKE_SYSTEM_NAME MATCHES "Linux")
        IF (NOT GCC_VERSION VERSION_LESS 4.2)
          SET (THEIA_CXX_FLAGS "${THEIA_CXX_FLAGS} -march=native -mtune=native")
        ENDIF (NOT GCC_VERSION VERSION_LESS 4.2)
      ENDIF (CMAKE_SYSTEM_NAME MATCHES "Linux")
      # Mac OS X
      IF (CMAKE_SYSTEM_NAME MATCHES "Darwin")
        SET (THEIA_CXX_FLAGS "${THEIA_CXX_FLAGS} -msse3")
        # Use of -fast only applicable for Apple's GCC
        # Assume this is being used if GCC version < 4.3 on OSX
        EXECUTE_PROCESS(COMMAND ${CMAKE_C_COMPILER}
          ARGS ${CMAKE_CXX_COMPILER_ARG1} -dumpversion
          OUTPUT_VARIABLE GCC_VERSION
          OUTPUT_STRIP_TRAILING_WHITESPACE)
        IF (GCC_VERSION VERSION_LESS 4.3)
          SET (THEIA_CXX_FLAGS "${THEIA_CXX_FLAGS} -fast")
        ENDIF (GCC_VERSION VERSION_LESS 4.3)
      ENDIF (CMAKE_SYSTEM_NAME MATCHES "Darwin")
    ENDIF (CMAKE_COMPILER_IS_GNUCXX)
    IF (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
      # Use of -O3 requires use of gold linker & LLVM-gold plugin, which might
      # well not be present / in use and without which files will compile, but
      # not link ('file not recognized') so explicitly check for support
      INCLUDE(CheckCXXCompilerFlag)
      CHECK_CXX_COMPILER_FLAG("-flto" HAVE_LTO_SUPPORT)
      IF (HAVE_LTO_SUPPORT)
        MESSAGE(STATUS "Enabling link-time optimization (-flto)")
        SET(THEIA_CXX_FLAGS "${THEIA_CXX_FLAGS} -flto")
      ELSE ()
        MESSAGE(STATUS "Compiler/linker does not support link-time optimization (-flto), disabling.")
      ENDIF (HAVE_LTO_SUPPORT)
    ENDIF ()
  ENDIF (CMAKE_BUILD_TYPE STREQUAL "Release")

  # Set c++ standard to c++11
  if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
      # Mac OS X before Mavericks uses libstdc++ by default but does not support
      # c++11. Force it to use libc++.
      IF (CMAKE_SYSTEM_NAME MATCHES "Darwin")
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
      ENDIF (CMAKE_SYSTEM_NAME MATCHES "Darwin")
  else()
    execute_process(COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
    if (GCC_VERSION VERSION_GREATER 4.7 OR GCC_VERSION VERSION_EQUAL 4.7)
      SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
    elseif(GCC_VERSION VERSION_GREATER 4.3 OR GCC_VERSION VERSION_EQUAL 4.3)
      message(WARNING "C++0x activated. If you get any errors update to a compiler which fully supports C++11")
      add_definitions("-std=gnu++0x")
      SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++0x")
    else ()
      message(FATAL_ERROR "C++11 needed. Therefore a gcc compiler with a version higher than 4.3 is needed.")
    endif()
  endif()

  SET (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${THEIA_CXX_FLAGS}")

  # After the tweaks for the compile settings, disable some warnings on MSVC.
  IF (MSVC)
    # Disable signed/unsigned int conversion warnings.
    ADD_DEFINITIONS("/wd4018")
    # Disable warning about using struct/class for the same symobl.
    ADD_DEFINITIONS("/wd4099")
    # Disable warning about the insecurity of using "std::copy".
    ADD_DEFINITIONS("/wd4996")
    # Disable performance warning about int-to-bool conversion.
    ADD_DEFINITIONS("/wd4800")
    # Disable performance warning about fopen insecurity.
    ADD_DEFINITIONS("/wd4996")
    # Disable warning about int64 to int32 conversion. Disabling
    # this warning may not be correct; needs investigation.
    # TODO(keir): Investigate these warnings in more detail.
    ADD_DEFINITIONS("/wd4244")
    # It's not possible to use STL types in DLL interfaces in a portable and
    # reliable way. However, that's what happens with Google Log and Google Flags
    # on Windows. MSVC gets upset about this and throws warnings that we can't do
    # much about. The real solution is to link static versions of Google Log and
    # Google Test, but that seems tricky on Windows. So, disable the warning.
    ADD_DEFINITIONS("/wd4251")

    # Google Flags doesn't have their DLL import/export stuff set up correctly,
    # which results in linker warnings. This is irrelevant for Theia, so ignore
    # the warnings.
    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /ignore:4049")

    # Update the C/CXX flags for MSVC to use either the static or shared
    # C-Run Time (CRT) library based on the user option: MSVC_USE_STATIC_CRT.
    LIST(APPEND C_CXX_FLAGS
      CMAKE_CXX_FLAGS
      CMAKE_CXX_FLAGS_DEBUG
      CMAKE_CXX_FLAGS_RELEASE
      CMAKE_CXX_FLAGS_MINSIZEREL
      CMAKE_CXX_FLAGS_RELWITHDEBINFO)

    FOREACH(FLAG_VAR ${C_CXX_FLAGS})
      IF (MSVC_USE_STATIC_CRT)
        # Use static CRT.
        IF (${FLAG_VAR} MATCHES "/MD")
          STRING(REGEX REPLACE "/MD" "/MT" ${FLAG_VAR} "${${FLAG_VAR}}")
        ENDIF (${FLAG_VAR} MATCHES "/MD")
      ELSE (MSVC_USE_STATIC_CRT)
        # Use shared, not static, CRT.
        IF (${FLAG_VAR} MATCHES "/MT")
          STRING(REGEX REPLACE "/MT" "/MD" ${FLAG_VAR} "${${FLAG_VAR}}")
        ENDIF (${FLAG_VAR} MATCHES "/MT")
      ENDIF (MSVC_USE_STATIC_CRT)
    ENDFOREACH()

    # Tuple sizes of 10 are used by Gtest.
    ADD_DEFINITIONS("-D_VARIADIC_MAX=10")
  ENDIF (MSVC)

  IF (UNIX)
    # GCC is not strict enough by default, so enable most of the warnings.
    SET(CMAKE_CXX_FLAGS
      "${CMAKE_CXX_FLAGS} -Werror=all -Werror=extra -Wno-unknown-pragmas -Wno-sign-compare -Wno-unused-parameter -Wno-missing-field-initializers")
  ENDIF (UNIX)

endmacro(OptimizeTheiaCompilerFlags)
