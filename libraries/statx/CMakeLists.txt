# Copyright (C) 2014  Victor Fragoso <vfragoso@cs.ucsb.edu>
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
#     * Neither the name of the University of California, Santa Barbara nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL VICTOR FRAGOSO BE LIABLE FOR ANY DIRECT,
# INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#

cmake_minimum_required(VERSION 2.8.0)
cmake_policy(VERSION 2.8)

if (COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif (COMMAND cmake_policy)

project (statx C CXX)
include_directories(${PROJECT_SOURCE_DIR})

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/libraries/statx/cmake/Modules/")

# C++11
# Compiler-specific C++11 activation.
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
    message(FATAL_ERROR "A compiler with C++11 support is required for Theia.")
  endif (COMPILER_HAS_CXX11_FLAG)
endif (NOT MSVC)
message("-- CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_INSTALL_PREFIX /usr/local/)

message("-- Default Install prefix: ${CMAKE_INSTALL_PREFIX}")

# Activate CERES
option(STATX_WITH_CERES "Enables GEV parameter estimation using Ceres" ON)

# Eigen
find_package(Eigen REQUIRED)
if (EIGEN_FOUND)
  message("-- Found Eigen version ${EIGEN_VERSION}: ${EIGEN_INCLUDE_DIRS}")
endif (EIGEN_FOUND)
include_directories(${EIGEN_INCLUDE_DIRS})

# Google Flags
find_package(Gflags REQUIRED)
if(GFLAGS_FOUND)
include_directories(${GFLAGS_INCLUDE_DIRS})
endif(GFLAGS_FOUND)

# Google Logging
find_package(Glog REQUIRED)
if(GLOG_FOUND)
include_directories(${GLOG_INCLUDE_DIRS})
endif(GLOG_FOUND)

# Google ceres-solver
if(STATX_WITH_CERES)
  # Check for CERES
  message("-- Check for Google ceres-solver")
  find_package(Ceres REQUIRED)
  if (CERES_FOUND)
    message("-- Found Ceres headers in: ${CERES_INCLUDE_DIRS}")
    include_directories(${CERES_INCLUDE_DIR})
  else (CERES_FOUND)
    message(FATAL_ERROR "Could not find Ceres. Please set CERES_INCLUDE_DIRS and CERES_LIBRARIES")
  endif(CERES_FOUND)

  # BLAS
  find_package(BLAS REQUIRED)

  # LAPACK
  find_package(LAPACK REQUIRED)
endif(STATX_WITH_CERES)

# Find and include optimo headers
message("-- Check for Optimo")
find_path(OPTIMO_INCLUDE_DIR
  NAMES optimo/solvers/bfgs.h
  PATHS ${CMAKE_SOURCE_DIR}/libraries/optimo /usr/include /usr/local/include)

if(OPTIMO_INCLUDE_DIR)
  message("-- Found optimo headers: ${OPTIMO_INCLUDE_DIR}")
else(OPTIMO_INCLUDE_DIR)
  message(FATAL_ERROR "Can't find optimo headers. Try passing -DOPTIMO_INCLUDE_DIR=...")
endif(OPTIMO_INCLUDE_DIR)
include_directories(${OPTIMO_INCLUDE_DIR})

# Polygamma functions from burkardt
set(XGAMMA_FUNCS ${CMAKE_CURRENT_SOURCE_DIR}/libraries/burkardt_spec_funcs)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/libraries/burkardt_spec_funcs)
set(XGAMMA_FILES ${XGAMMA_FUNCS}/asa121.cpp ${XGAMMA_FUNCS}/asa103.cpp)

# Setting CXX FLAGS appropriately. The code below was inspired from
# Google CERES and modified for this library.
set (CMAKE_BUILD_TYPE Release)
set (STATX_CXX_FLAGS)

if (CMAKE_BUILD_TYPE STREQUAL "Release")
  if (CMAKE_COMPILER_IS_GNUCXX)
    # Linux
    if (CMAKE_SYSTEM_NAME MATCHES "Linux")
      if (NOT GCC_VERSION VERSION_LESS 4.2)
        set (STATX_CXX_FLAGS
          "${STATX_CXX_FLAGS} -march=native -mtune=native -msse2 -msse3 -msse4")
      endif (NOT GCC_VERSION VERSION_LESS 4.2)
    endif (CMAKE_SYSTEM_NAME MATCHES "Linux")
    # Mac OS X
    if (CMAKE_SYSTEM_NAME MATCHES "Darwin")
      set (STATX_CXX_FLAGS "${STATX_CXX_FLAGS} -msse3 -msse4 -msse3 -msse2")
      execute_process(COMMAND ${CMAKE_C_COMPILER}
        ARGS ${CMAKE_CXX_COMPILER_ARG1} -dumpversion
        OUTPUT_VARIABLE GCC_VERSION
        OUTPUT_STRIP_TRAILING_WHITESPACE)
      if (GCC_VERSION VERSION_LESS 4.3)
        set (STATX_CXX_FLAGS "${STATX_CXX_FLAGS} -fast")
      endif (GCC_VERSION VERSION_LESS 4.3)
    endif (CMAKE_SYSTEM_NAME MATCHES "Darwin")
  endif (CMAKE_COMPILER_IS_GNUCXX)
  if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    set (STATX_CXX_FLAGS "-Ofast -ffast-math -fvectorize -funroll-loops")
  endif ()
endif (CMAKE_BUILD_TYPE STREQUAL "Release")

set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${STATX_CXX_FLAGS} -DNDEBUG")

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
  if (CMAKE_COMPILER_IS_GNUCXX)
    # Linux
    if (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
      set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g3 -O0")
    endif (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    # Mac OS X
    if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
      set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g3 -O0")
    endif (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  endif (CMAKE_COMPILER_IS_GNUCXX)
  # CLANG
  if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_DEBUG} -g3 -O0")
  endif ()
endif (CMAKE_BUILD_TYPE STREQUAL "Debug")

# Main Library
# Note: Gumbel parameter estimator was removed due to compatibility issues with Eigen 3.2.3.
set(STATX_DISTRIBUTIONS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/statx/distributions)
set(EVD_PREFIX_DIR ${STATX_DISTRIBUTIONS_DIR}/evd)
set(UTILS_PREFIX_DIR statx/utils)
set(statx_src ${UTILS_PREFIX_DIR}/ecdf.cc ${XGAMMA_FILES}
  ${STATX_DISTRIBUTIONS_DIR}/gamma.cc
  ${EVD_PREFIX_DIR}/gev_mle.cc ${EVD_PREFIX_DIR}/gev.cc
  ${EVD_PREFIX_DIR}/gpd_mle.cc ${EVD_PREFIX_DIR}/gpd.cc)
set(default_link_libraries ${GLOG_LIBRARIES})

# Check if CERES was activated
if(STATX_WITH_CERES)
set(statx_src ${statx_src} ${EVD_PREFIX_DIR}/gev_ceres.cc ${EVD_PREFIX_DIR}/gpd_ceres.cc)
add_definitions(-DSTATX_WITH_CERES)
set(default_link_libraries ${default_link_libraries} ${CERES_LIBRARIES}
  ${SUITESPARSE_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
endif(STATX_WITH_CERES)
add_library(statx STATIC ${statx_src})

target_link_libraries(statx ${default_link_libraries})

# Adding subdirectory
add_subdirectory(statx)

# Exporting target for Theia
install(TARGETS statx
  EXPORT  TheiaExport
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib${LIB_SUFFIX}
  ARCHIVE DESTINATION lib${LIB_SUFFIX})
