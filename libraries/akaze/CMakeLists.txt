cmake_minimum_required(VERSION 2.8)
project( akaze )
set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)

OPTION(OPENMP "Enable multithreading (requires OpenMP)" OFF)

# Where to output binaries and libraries
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/lib")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin")

FIND_PACKAGE(Eigen REQUIRED)
INCLUDE_DIRECTORIES(${EIGEN_INCLUDE_DIRS})
IF (EIGEN_FOUND)
  MESSAGE("-- Found Eigen version ${EIGEN_VERSION}: ${EIGEN_INCLUDE_DIRS}")
ENDIF (EIGEN_FOUND)

# CImg is needed for reading/displaying the images.
ADD_SUBDIRECTORY(cimg)
#INCLUDE_DIRECTORIES(${cimg_INCLUDE_DIR} ${timer_INCLUDE_DIR} ./)
INCLUDE_DIRECTORIES(./)

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

IF (OPENMP)
  # Clang does not (yet) support OpenMP.
  IF (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    # Retain the help string associated with the OPENMP option
    # when updating it to disable use of OPENMP.
    GET_PROPERTY(HELP_STRING CACHE OPENMP PROPERTY HELPSTRING)
    SET(OPENMP OFF CACHE BOOL "${HELP_STRING}" FORCE)
    MESSAGE("-- Compiler is Clang, disabling OpenMP.")
  ELSE (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    # Find quietly s/t as we can continue without OpenMP if it is not found.
    FIND_PACKAGE(OpenMP QUIET)
    IF (OPENMP_FOUND)
      MESSAGE("-- Building with OpenMP.")
      ADD_DEFINITIONS(-DAKAZE_USE_OPENMP)
      SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
      IF (UNIX)
        # At least on Linux, we need pthreads to be enabled for mutex to
        # compile.  This may not work on Windows or Android.
        FIND_PACKAGE(Threads REQUIRED)
        SET(STATIC_LIBRARY_FLAGS
          "${STATIC_LIBRARY_FLAGS} ${CMAKE_THREAD_LIBS_INIT}")
        SET(CMAKE_SHARED_LINKER_FLAGS
          "${CMAKE_SHARED_LINKER_FLAGS} ${CMAKE_THREAD_LIBS_INIT}")
      ENDIF (UNIX)
    ELSE (OPENMP_FOUND)
      MESSAGE("-- Failed to find OpenMP, disabling.")
      # Retain the help string associated with the OPENMP option
      # when updating it to disable use of OPENMP.
      GET_PROPERTY(HELP_STRING CACHE OPENMP PROPERTY HELPSTRING)
      SET(OPENMP OFF CACHE BOOL "${HELP_STRING}" FORCE)
    ENDIF (OPENMP_FOUND)
  ENDIF (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
ELSE (OPENMP)
  MESSAGE("-- Building without OpenMP.")
ENDIF (OPENMP)

Include("${PROJECT_SOURCE_DIR}/cmake/OptimizeCompilerFlags.cmake")
OptimizeCompilerFlags()

# Support for Debug mode.
IF (CMAKE_BUILD_TYPE STREQUAL "Debug")
  MESSAGE("========== Warning: Compiling with Debug mode ==========")
ENDIF (CMAKE_BUILD_TYPE STREQUAL "Debug")

mark_as_advanced(
  Gperftools_ROOT_DIR
  GPERFTOOLS_TCMALLOC
  GPERFTOOLS_PROFILER
  GPERFTOOLS_TCMALLOC_AND_PROFILER
  GPERFTOOLS_LIBRARIES
  GPERFTOOLS_INCLUDE_DIR)

SET(AKAZE_SRCS
    src/AKAZE.cpp
    src/convolution.cpp
    src/fed.cpp
    src/nldiffusion_functions.cpp
    src/utils.cpp
    timer/timer.cpp)
ADD_LIBRARY(akaze ${AKAZE_SRCS})
TARGET_LINK_LIBRARIES(akaze ${EIGEN_LIBRARIES} ${CIMG_LIBRARIES})

# Demo program extracting akaze features.
ADD_EXECUTABLE(akaze_features akaze_features.cpp)
TARGET_LINK_LIBRARIES(akaze_features akaze ${CIMG_LIBRARIES})

# Demo program matching akaze features.
ADD_EXECUTABLE(akaze_match akaze_match.cpp)
TARGET_LINK_LIBRARIES(akaze_match akaze ${CIMG_LIBRARIES})

install(TARGETS akaze
  EXPORT  TheiaExport
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib${LIB_SUFFIX}
  ARCHIVE DESTINATION lib${LIB_SUFFIX})
