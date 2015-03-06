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

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.3)
IF (COMMAND cmake_policy)
  CMAKE_POLICY(SET CMP0003 NEW)
ENDIF (COMMAND cmake_policy)

# Disable warnings about rpath on MacOS X for now.
# TODO(cmsweeney): Figure out the proper way to handle this!
IF(POLICY CMP0042)
  CMAKE_POLICY(SET CMP0042 NEW) # CMake 3.0
ENDIF()

PROJECT(THEIA C CXX)
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR})

# Set up the git hook to make Gerrit Change-Id: lines in commit messages.
SET (LOCAL_GIT_DIRECTORY)
IF (EXISTS ${CMAKE_SOURCE_DIR}/.git)
  # .git directory can be found on Unix based system, or on Windows with
  # Git Bash (shipped with msysgit)
  SET (LOCAL_GIT_DIRECTORY ${CMAKE_SOURCE_DIR}/.git)
ELSE (EXISTS ${CMAKE_SOURCE_DIR}/.git)
ENDIF (EXISTS ${CMAKE_SOURCE_DIR}/.git)

IF (EXISTS ${LOCAL_GIT_DIRECTORY})
  IF (NOT EXISTS ${LOCAL_GIT_DIRECTORY}/hooks/commit-msg)
    # Download the hook only if it is not already present. We borrow the hook
    # hosted by Ceres because it is not specific to the project.
    FILE(DOWNLOAD http://cs.ucsb.edu/~cmsweeney/theia/commit-msg
         ${CMAKE_BINARY_DIR}/commit-msg)

    # Make the downloaded file executable, since it is not by default.
    FILE(COPY ${CMAKE_BINARY_DIR}/commit-msg
         DESTINATION ${LOCAL_GIT_DIRECTORY}/hooks/
         FILE_PERMISSIONS
           OWNER_READ OWNER_WRITE OWNER_EXECUTE
           GROUP_READ GROUP_WRITE GROUP_EXECUTE
           WORLD_READ WORLD_EXECUTE)
  ENDIF (NOT EXISTS ${LOCAL_GIT_DIRECTORY}/hooks/commit-msg)
ENDIF (EXISTS ${LOCAL_GIT_DIRECTORY})

# ==============================================================================
# Additional cmake find modules
# ==============================================================================
LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

OPTION(BUILD_TESTING "Enable testing" ON)
OPTION(BUILD_DOCUMENTATION "Build html User's Guide" OFF)

ENABLE_TESTING()
ADD_DEFINITIONS(-DGTEST_USE_OWN_TR1_TUPLE=1)
SET(GTEST_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/libraries/gtest")
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

SET(CMAKE_RELEASE_POSTFIX "")
SET(CMAKE_DEBUG_POSTFIX "-debug")

SET(THEIA_VERSION_MAJOR 0)
SET(THEIA_VERSION_MINOR 1)
SET(THEIA_VERSION_PATCH 0)
SET(THEIA_VERSION
    ${THEIA_VERSION_MAJOR}.${THEIA_VERSION_MINOR}.${THEIA_VERSION_PATCH})
SET(THEIA_ABI_VERSION 0.1.0)

# THEIA data directory
ADD_DEFINITIONS(-DTHEIA_DATA_DIR="${CMAKE_SOURCE_DIR}/data")

# Eigen
SET(MIN_EIGEN_VERSION 3.2.0)
FIND_PACKAGE(Eigen ${MIN_EIGEN_VERSION} REQUIRED)
IF (EIGEN_FOUND)
  MESSAGE("-- Found Eigen version ${EIGEN_VERSION}: ${EIGEN_INCLUDE_DIRS}")
ENDIF (EIGEN_FOUND)

# Use a larger inlining threshold for Clang, since it hobbles Eigen,
# resulting in an unreasonably slow version of the blas routines. The
# -Qunused-arguments is needed because CMake passes the inline
# threshold to the linker and clang complains about it and dies.
IF (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  SET(CMAKE_CXX_FLAGS
    "${CMAKE_CXX_FLAGS} -Qunused-arguments -mllvm -inline-threshold=600")
  # Older versions of Clang (<= 2.9) do not support the 'return-type-c-linkage'
  # option, so check for its presence before adding it to the default flags set.
  INCLUDE(CheckCXXCompilerFlag)
  CHECK_CXX_COMPILER_FLAG("-Wno-return-type-c-linkage"
    HAVE_RETURN_TYPE_C_LINKAGE)
  IF (HAVE_RETURN_TYPE_C_LINKAGE)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-return-type-c-linkage")
  ENDIF(HAVE_RETURN_TYPE_C_LINKAGE)
ENDIF ()

# GFlags. The namespace patch is borrow from Ceres Solver (see license in
# FindGflags.cmake)
FIND_PACKAGE(Gflags REQUIRED)
IF (GFLAGS_FOUND)
  MESSAGE("-- Found Google Flags: ${GFLAGS_INCLUDE_DIRS} in namespace: ${GFLAGS_NAMESPACE}")
  ADD_DEFINITIONS(-DTHEIA_GFLAGS_NAMESPACE=${GFLAGS_NAMESPACE})
ELSE (GLOG_FOUND)
  MESSAGE(FATAL_ERROR "Can't find Google FLAGS.")
ENDIF (GFLAGS_FOUND)

# Google Logging
MESSAGE("-- Check for Google Log")
FIND_PACKAGE(Glog REQUIRED)
IF (GLOG_FOUND)
  MESSAGE("-- Found Google Logging: ${GLOG_INCLUDE_DIRS}")
ELSE (GLOG_FOUND)
  MESSAGE(FATAL_ERROR "Can't find Google Logging. Please set GLOG_INCLUDE_DIR & "
    "GLOG_LIBRARY")
ENDIF (GLOG_FOUND)

# Ceres
MESSAGE("-- Check for Ceres")
FIND_PACKAGE(Ceres REQUIRED)
IF (CERES_FOUND)
  MESSAGE("-- Found Ceres: ${CERES_INCLUDE_DIRS}")
ELSE (CERES_FOUND)
  MESSAGE(FATAL_ERROR "Can't find Ceres. Please set CERES_INCLUDE_DIR & "
    "CERES_LIBRARY")
ENDIF (CERES_FOUND)

INCLUDE_DIRECTORIES(
  include
  src
  ${CERES_INCLUDE_DIR}
  ${GFLAGS_INCLUDE_DIRS}
  ${GLOG_INCLUDE_DIRS})

# NOTE: This fix came from Ceres solver with the following comment:
#
# Eigen SparseQR generates various compiler warnings related to unused and
# uninitialised local variables, which prevents Ceres compilation as we use
# -Werror.  To avoid having to individually suppress these warnings around
# the #include statments for Eigen headers across all GCC/Clang versions, we
# tell CMake to treat Eigen headers as system headers.  This results in all
# compiler warnings from them being suppressed.
#
# Note that this is *not* propagated to clients, ie CERES_INCLUDE_DIRS
# used by clients after find_package(Ceres) does not identify Eigen as
# as system headers.
INCLUDE_DIRECTORIES(SYSTEM ${EIGEN_INCLUDE_DIRS})

# build Third party libraries included in distro.
ADD_SUBDIRECTORY(libraries)
INCLUDE_DIRECTORIES(
  libraries
  ${agast_SOURCE_DIR}
  ${cimg_SOURCE_DIR}
  ${cimg_INCLUDE_DIR}
  ${easyexif_INCLUDE_DIR}
  ${gtest_SOURCE_DIR}/include
  ${gtest_SOURCE_DIR}
  ${optimo_SOURCE_DIR}
  ${statx_SOURCE_DIR}
  ${stlplus3_SOURCE_DIR}
  ${vlfeat_SOURCE_DIR})

LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
INCLUDE(OptimizeTheiaCompilerFlags)
OptimizeTheiaCompilerFlags()

ADD_SUBDIRECTORY(src/theia)
ADD_SUBDIRECTORY(applications)

IF (BUILD_DOCUMENTATION)
  MESSAGE("-- Documentation building is enabled")

  # Generate the User's Guide (html).
  # The corresponding target is UserGuide, but is included in ALL.
  ADD_SUBDIRECTORY(docs)
ENDIF (BUILD_DOCUMENTATION)

# ---------------- Install commands ----------------- #
# Install the .h files
FILE(GLOB THEIA_HDRS ${CMAKE_SOURCE_DIR}/include/theia/*.h)
INSTALL(FILES ${THEIA_HDRS} DESTINATION include/theia)

FILE(GLOB_RECURSE THEIA_INTERNAL_HDRS ${CMAKE_SOURCE_DIR}/src/theia/*.h)
INSTALL(DIRECTORY src/theia/ DESTINATION include/theia FILES_MATCHING PATTERN "*.h")

INSTALL(DIRECTORY libraries/ DESTINATION include/theia/libraries FILES_MATCHING PATTERN "*.h*")

# Add an uninstall target to remove all installed files.
CONFIGURE_FILE("${CMAKE_SOURCE_DIR}/cmake/uninstall.cmake.in"
               "${CMAKE_BINARY_DIR}/cmake/uninstall.cmake"
               IMMEDIATE @ONLY)

ADD_CUSTOM_TARGET(uninstall
                  COMMAND ${CMAKE_COMMAND} -P ${CMAKE_BINARY_DIR}/cmake/uninstall.cmake)

# Set up install directories. INCLUDE_INSTALL_DIR, LIB_INSTALL_DIR and
# CMAKECONFIG_INSTALL_DIR must not be absolute paths.
IF (WIN32)
  SET(INCLUDE_INSTALL_DIR Include)
  SET(LIB_INSTALL_DIR Lib)
  SET(CMAKECONFIG_INSTALL_DIR CMake)
  SET(RELATIVE_CMAKECONFIG_INSTALL_DIR CMake)
ELSE ()
  SET(INCLUDE_INSTALL_DIR include)
  SET(LIB_INSTALL_DIR lib)
  SET(CMAKECONFIG_INSTALL_DIR share/Theia)
  SET(RELATIVE_CMAKECONFIG_INSTALL_DIR share/Theia)
ENDIF ()

# This "exports" all targets which have been put into the export set
# "TheiaExport". This means that CMake generates a file with the given
# filename, which can later on be loaded by projects using this package.
# This file contains ADD_LIBRARY(bar IMPORTED) statements for each target
# in the export set, so when loaded later on CMake will create "imported"
# library targets from these, which can be used in many ways in the same way
# as a normal library target created via a normal ADD_LIBRARY().
INSTALL(EXPORT TheiaExport
        DESTINATION ${RELATIVE_CMAKECONFIG_INSTALL_DIR} FILE TheiaTargets.cmake)

# Figure out the relative path from the installed Config.cmake file to the
# install prefix (which may be at runtime different from the chosen
# CMAKE_INSTALL_PREFIX if under Windows the package was installed anywhere)
# This relative path will be configured into the TheiaConfig.cmake.
FILE(RELATIVE_PATH INSTALL_ROOT_REL_CONFIG_INSTALL_DIR
     ${CMAKE_INSTALL_PREFIX}/${CMAKECONFIG_INSTALL_DIR} ${CMAKE_INSTALL_PREFIX})

# Create a TheiaConfig.cmake file. <name>Config.cmake files are searched by
# FIND_PACKAGE() automatically. We configure that file so that we can put any
# information we want in it, e.g. version numbers, include directories, etc.
CONFIGURE_FILE("${CMAKE_SOURCE_DIR}/cmake/TheiaConfig.cmake.in"
               "${CMAKE_CURRENT_BINARY_DIR}/TheiaConfig.cmake" @ONLY)

# Additionally, when CMake has found a TheiaConfig.cmake, it can check for a
# TheiaConfigVersion.cmake in the same directory when figuring out the version
# of the package when a version has been specified in the FIND_PACKAGE() call,
# e.g. FIND_PACKAGE(Theia [0.5.2] REQUIRED). The version argument is optional.
CONFIGURE_FILE("${CMAKE_SOURCE_DIR}/cmake/TheiaConfigVersion.cmake.in"
               "${CMAKE_CURRENT_BINARY_DIR}/TheiaConfigVersion.cmake" @ONLY)

# Install these files into the same directory as the generated exports-file,
# we include the FindPackage scripts for libraries whose headers are included
# in the public API of Theia and should thus be present in THEIA_INCLUDE_DIRS.
INSTALL(FILES "${CMAKE_CURRENT_BINARY_DIR}/TheiaConfig.cmake"
              "${CMAKE_CURRENT_BINARY_DIR}/TheiaConfigVersion.cmake"
              "${CMAKE_SOURCE_DIR}/cmake/FindEigen.cmake"
              "${CMAKE_SOURCE_DIR}/cmake/FindGlog.cmake"
              "${CMAKE_SOURCE_DIR}/cmake/FindGflags.cmake"
              DESTINATION ${CMAKECONFIG_INSTALL_DIR})
