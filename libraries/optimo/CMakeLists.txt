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

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.0)
CMAKE_POLICY(VERSION 2.8)

IF (COMMAND cmake_policy)
  CMAKE_POLICY(SET CMP0003 NEW)
ENDIF (COMMAND cmake_policy)

PROJECT(optimo C CXX)
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR})

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# C++11
# Compiler-specific C++11 activation.
if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
    execute_process(
        COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
    if (NOT (GCC_VERSION VERSION_GREATER 4.7 OR GCC_VERSION VERSION_EQUAL 4.7))
        message(FATAL_ERROR "${PROJECT_NAME} requires g++ 4.7 or greater.")
    endif ()
    set(CMAKE_CXX_FLAGS "-std=c++11")
elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    set(CMAKE_CXX_FLAGS "-std=c++11 -stdlib=libc++")
else ()
    message(FATAL_ERROR "Your C++ compiler does not support C++11.")
endif ()
# Showing the C++ compiler flags
MESSAGE("-- CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Add GTest as a target.
OPTION(OPTIMO_BUILD_TESTING "Enable testing" OFF)
IF (${OPTIMO_BUILD_TESTING})
  ENABLE_TESTING()
  ADD_SUBDIRECTORY(libraries/gtest)
  INCLUDE_DIRECTORIES(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})
ENDIF (${OPTIMO_BUILD_TESTING})

# Locations to search for Eigen
# Eigen
FIND_PACKAGE(Eigen REQUIRED)
IF (EIGEN_FOUND)
  MESSAGE("-- Found Eigen version ${EIGEN_VERSION}: ${EIGEN_INCLUDE_DIRS}")
ENDIF (EIGEN_FOUND)
INCLUDE_DIRECTORIES(${EIGEN_INCLUDE_DIRS})

# Setting CXX FLAGS appropriately. The code below was inspired from
# Google CERES and modified for this library.
SET (CMAKE_BUILD_TYPE Release)
SET (OPTIMO_CXX_FLAGS)

IF (CMAKE_BUILD_TYPE STREQUAL "Release")
  IF (CMAKE_COMPILER_IS_GNUCXX)
    # Linux
    IF (CMAKE_SYSTEM_NAME MATCHES "Linux")
      IF (NOT GCC_VERSION VERSION_LESS 4.2)
        SET (OPTIMO_CXX_FLAGS "${OPTIMO_CXX_FLAGS} -march=native -mtune=native -msse2 -msse3 -msse4")
      ENDIF (NOT GCC_VERSION VERSION_LESS 4.2)
    ENDIF (CMAKE_SYSTEM_NAME MATCHES "Linux")
    # Mac OS X
    IF (CMAKE_SYSTEM_NAME MATCHES "Darwin")
      SET (OPTIMO_CXX_FLAGS "${OPTIMO_CXX_FLAGS} -msse3 -msse4 -msse3 -msse2")
      EXECUTE_PROCESS(COMMAND ${CMAKE_C_COMPILER}
        ARGS ${CMAKE_CXX_COMPILER_ARG1} -dumpversion
        OUTPUT_VARIABLE GCC_VERSION
        OUTPUT_STRIP_TRAILING_WHITESPACE)
      IF (GCC_VERSION VERSION_LESS 4.3)
        SET (OPTIMO_CXX_FLAGS "${OPTIMO_CXX_FLAGS} -fast")
      ENDIF (GCC_VERSION VERSION_LESS 4.3)
    ENDIF (CMAKE_SYSTEM_NAME MATCHES "Darwin")
  ENDIF (CMAKE_COMPILER_IS_GNUCXX)
  IF (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    SET (OPTIMO_CXX_FLAGS "-Ofast -ffast-math -fvectorize -funroll-loops")
  ENDIF ()
ENDIF (CMAKE_BUILD_TYPE STREQUAL "Release")

SET (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${OPTIMO_CXX_FLAGS} -DNDEBUG")

IF (CMAKE_BUILD_TYPE STREQUAL "Debug")
  IF (CMAKE_COMPILER_IS_GNUCXX)
    # Linux
    IF (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
      SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g3 -O0")
    ENDIF (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    # Mac OS X
    IF (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
      SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g3 -O0")
    ENDIF (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  ENDIF (CMAKE_COMPILER_IS_GNUCXX)
  # CLANG
  IF (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_DEBUG} -g3 -O0")
  ENDIF ()
ENDIF (CMAKE_BUILD_TYPE STREQUAL "Debug")

IF (OPTIMO_BUILD_TESTING)
  # Google Flags
  FIND_PACKAGE(Gflags REQUIRED)
  IF(GFLAGS_FOUND)
    INCLUDE_DIRECTORIES(${GFLAGS_INCLUDE_DIRS})
  ENDIF(GFLAGS_FOUND)

  # Google Logging
  FIND_PACKAGE(Glog REQUIRED)
  IF(GLOG_FOUND)
    INCLUDE_DIRECTORIES(${GLOG_INCLUDE_DIRS})
  ENDIF(GLOG_FOUND)

  # This macro was taken from Chris Sweeney's libtheia library.
  # Macro for tests
  #   Ex: Your test is written in myclass_test.cc  and uses a class in myclass.cc
  #       which has been compiled as target: myclass then calling
  #       GTEST(myclass_test myclass) will create a test executable myclass_test
  #       by automatically compiling myclass_test.cc and linking the mytest target
  ADD_LIBRARY(test_main test/test_main.cc)
  TARGET_LINK_LIBRARIES(test_main gtest ${GLOG_LIBRARIES} ${GFLAGS_LIBRARIES} ${blas_LIBRARIES}
    ${SUITESPARSE_LIBRARIES})
  MACRO (GTEST NAME)
    ADD_EXECUTABLE(${NAME} ${NAME}.cc)
    TARGET_LINK_LIBRARIES(${NAME} test_main gtest ${ARGN} ${GFLAGS_LIBRARIES} ${GLOG_LIBRARIES}
      ${blas_LIBRARIES} ${SUITESPARSE_LIBRARIES})
    ADD_TEST(NAME ${NAME}
      COMMAND ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${NAME})
  ENDMACRO (GTEST)
ENDIF(OPTIMO_BUILD_TESTING)

ADD_SUBDIRECTORY(optimo)

# Installation of Headers
INSTALL(DIRECTORY ${PROJECT_SOURCE_DIR}/optimo
  DESTINATION ${CMAKE_INSTALL_PREFIX}
  FILES_MATCHING PATTERN "*.h"
  PATTERN "CMakeFiles" EXCLUDE)