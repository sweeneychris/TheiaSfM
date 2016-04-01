# CMakeFiles.txt for CImg
#
# Haz-Edine Assemlal <assemlal@gmail.com>

cmake_minimum_required(VERSION 2.8.0)

set(CIMG_PACKAGE_VERSION_MAJOR 1)
set(CIMG_PACKAGE_VERSION_MINOR 5)
set(CIMG_PACKAGE_VERSION_PATCH 0)

# register local modules
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" "${CMAKE_SOURCE_DIR}/cmake-modules")


project(cimg)

# Packaging defines
set(CPACK_GENERATOR TGZ)
set(CPACK_PACKAGE_VERSION_MAJOR ${CIMG_PACKAGE_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${CIMG_PACKAGE_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${CIMG_PACKAGE_VERSION_PATCH})
include(CPack)

set(PACKAGE "CImg")
set(PACKAGE_BUGREPORT "david.tschumperle@greyc.ensicaen.fr")
set(PACKAGE_NAME "CImg")
set(PACKAGE_VERSION "${CIMG_PACKAGE_VERSION_MAJOR}.${CIMG_PACKAGE_VERSION_MINOR}.${CIMG_PACKAGE_VERSION_PATCH}")
set(PACKAGE_STRING "${PACKAGE_NAME} ${PACKAGE_VERSION}")


option(CIMG_WITH_LAPACK	  "Use LAPACK library" ON)

# Try to find all image libraries. Throw an error if none can be found.
find_package(PNG)
find_package(JPEG)
find_package(TIFF)

if (NOT PNG_FOUND AND NOT JPEG_FOUND AND NOT TIFF_FOUND)
  message(FATAL_ERROR "Cannot find libpng, libjpeg, or libtiff. At least one must be installed in order to use CImg for loading images.")
endif (NOT PNG_FOUND AND NOT JPEG_FOUND AND NOT TIFF_FOUND)

# ImageMagick is required for reading jpeg and png files.
message("-- Check for ImageMagick")
find_package(ImageMagick COMPONENTS convert mogrify REQUIRED)
if (ImageMagick_FOUND)
  message("-- Found ImageMagick: ${ImageMagick_INCLUDE_DIRS}")
  include_directories( ${ImageMagick_convert_INCLUDE_DIRS} ${ImageMagick_mogrify_INCLUDE_DIRS} )
  list( APPEND DEPENDENCIES_LIBRARIES ${ImageMagick_convert_LIBRARIES} ${ImageMagick_mogrify_LIBRARIES} )
else (ImageMagick_FOUND)
  message(FATAL_ERROR "Can't find ImageMagick. Please set ImageMagick_INCLUDE_DIR & "
    "ImageMagick_LIBRARY")
endif (ImageMagick_FOUND)

if(PNG_FOUND)
  add_definitions(-Dcimg_use_png ${PNG_DEFINITIONS})
endif(PNG_FOUND)

if(JPEG_FOUND)
  add_definitions( -Dcimg_use_jpeg )
endif(JPEG_FOUND)

if(TIFF_FOUND)
  add_definitions( -Dcimg_use_tiff )
endif(TIFF_FOUND)

if(CIMG_WITH_LAPACK)
  find_package( LAPACK QUIET )
  if (LAPACK_FOUND)
    add_definitions( -Dcimg_use_lapack )
    include_directories( ${LAPACK_INCLUDE_DIRS} )
    list( APPEND DEPENDENCIES_LIBRARIES ${LAPACK_LIBRARIES} )
  else (LAPACK_FOUND)
    message("-- Failed to find Lapack, disabling.")
  endif (LAPACK_FOUND)
endif(CIMG_WITH_LAPACK)

if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  if (CMAKE_SYSTEM_NAME MATCHES "Darwin")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-shift-negative-value")
  endif (CMAKE_SYSTEM_NAME MATCHES "Darwin")
endif (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")

# aliases
set(VERSION "${PACKAGE_VERSION}")
include_directories(  ${CMAKE_CURRENT_SOURCE_DIR} )
set(CIMG_LIBRARIES ${DEPENDENCIES_LIBRARIES} CACHE INTERNAL "include cimg libs")
