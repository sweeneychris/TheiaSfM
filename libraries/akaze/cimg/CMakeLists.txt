# CMakeFiles.txt for CImg
#
# Haz-Edine Assemlal <assemlal@gmail.com>

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.0)

set(CIMG_PACKAGE_VERSION_MAJOR 1)
set(CIMG_PACKAGE_VERSION_MINOR 5)
set(CIMG_PACKAGE_VERSION_PATCH 0)

# register local modules
SET(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" "${CMAKE_SOURCE_DIR}/cmake-modules")

PROJECT(cimg)

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


option(BUILD_GIMP	  "Build GIMP plug-ins" OFF)
option(BUILD_GMIC	  "Build GMIC" OFF)
option(WITH_VT100	  "Enable color output messages " ON)
option(WITH_CIMG_DISPLAY  "Use CImg display capabilities" OFF)
option(WITH_OPENMP	  "Use OpenMP parallel core instructions" ON)
option(WITH_XRANDR	  "Use XRandR library" OFF)
option(WITH_PNG		  "Use PNG library" OFF)
option(WITH_JPEG	  "Use JPEG library" ON)
option(WITH_TIFF	  "Use TIFF library" OFF)
option(WITH_ZLIB	  "Use Zlib library" OFF)
option(WITH_MAGICK	  "Use Magick library" OFF)
option(WITH_LAPACK	  "Use LAPACK library" OFF)

if(WITH_VT100)
  add_definitions(-Dcimg_use_vt100)
endif(WITH_VT100)

if( WITH_CIMG_DISPLAY )
  if(UNIX OR APPLE)
    if( CMAKE_SYSTEM_NAME MATCHES "(Solaris|SunOS)" )
       list( APPEND DEPENDENCIES_LIBRARIES "-R$(X11PATH)/lib -lrt -lnsl -lsocket" )
     else()
       find_package( X11 REQUIRED )
       include_directories( ${X11_INCLUDE_DIR} )
       list( APPEND DEPENDENCIES_LIBRARIES ${X11_LIBRARIES} )
     endif()
  elseif(WIN32)
    set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mwindows" )
    list( APPEND DEPENDENCIES_LIBRARIES "-lgdi32" )
  endif()
else(WITH_CIMG_DISPLAY)
  add_definitions(-Dcimg_display=0)
endif(WITH_CIMG_DISPLAY)

IF (WITH_OPENMP)
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
      ADD_DEFINITIONS(-Dcimg_use_openmp)
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
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
ELSE (WITH_OPENMP)
  MESSAGE("-- Building without OpenMP.")
ENDIF (WITH_OPENMP)

if(WITH_XRANDR)
  if(NOT X11_Xrandr_INCLUDE_PATH)
    message(FATAL_ERROR "Failed to find X11 XrandR which is required")
  endif(NOT X11_Xrandr_INCLUDE_PATH)
endif(WITH_XRANDR)

if(WITH_PNG)
  find_package( PNG REQUIRED )
  add_definitions(-Dcimg_use_png ${PNG_DEFINITIONS})
  include_directories( ${PNG_INCLUDE_DIRS} )
  list( APPEND DEPENDENCIES_LIBRARIES ${PNG_LIBRARIES} )
endif(WITH_PNG)

if(WITH_JPEG)
  find_package( JPEG REQUIRED )
  add_definitions( -Dcimg_use_jpeg )
  include_directories( ${JPEG_INCLUDE_DIR} )
  list( APPEND DEPENDENCIES_LIBRARIES ${JPEG_LIBRARIES} )
endif(WITH_JPEG)

if(WITH_TIFF)
  find_package( TIFF REQUIRED )
  add_definitions( -Dcimg_use_tiff )
  include_directories( ${TIFF_INCLUDE_DIR} )
  list( APPEND DEPENDENCIES_LIBRARIES ${TIFF_LIBRARIES} )
endif(WITH_TIFF)

if(WITH_ZLIB)
  find_package( ZLIB REQUIRED )
  add_definitions( -Dcimg_use_zlib )
  include_directories( ${ZLIB_INCLUDE_DIR} )
  list( APPEND DEPENDENCIES_LIBRARIES ${ZLIB_LIBRARIES} )
endif(WITH_ZLIB)

if(WITH_MAGICK)
  find_package( ImageMagick QUIET COMPONENTS Magick++ convert)
  IF (ImageMagick_FOUND)
    add_definitions( -Dcimg_use_magick )
    include_directories( ${ImageMagick_INCLUDE_DIRS} )
    list( APPEND DEPENDENCIES_LIBRARIES ${ImageMagick_LIBRARIES} )
  ELSE (ImageMagick_FOUND)
    MESSAGE("-- Failed to find ImageMagick, disabling.")
  ENDIF(ImageMagick_FOUND)
endif(WITH_MAGICK)

if(WITH_LAPACK)
  find_package( LAPACK QUIET )
  IF (LAPACK_FOUND)
    add_definitions( -Dcimg_use_lapack )
    include_directories( ${LAPACK_INCLUDE_DIRS} )
    list( APPEND DEPENDENCIES_LIBRARIES ${LAPACK_LIBRARIES} )
  ELSE (LAPACK_FOUND)
    MESSAGE("-- Failed to find Lapack, disabling.")
  ENDIF (LAPACK_FOUND)
endif(WITH_LAPACK)

if(WITH_BOARD)
  find_package( BOARD  )
  IF (BOARD_FOUND)
    add_definitions( -Dcimg_use_board )
    include_directories( ${BOARD_INCLUDE_DIR} )
    list( APPEND DEPENDENCIES_LIBRARIES ${BOARD_LIBRARIES} )
  ELSE (BOARD_FOUND)
    MESSAGE("-- Failed to find Board library, disabling.")
  ENDIF (BOARD_FOUND)
endif(WITH_BOARD)

# aliases
set(VERSION "${PACKAGE_VERSION}")
include_directories(  ${CMAKE_CURRENT_SOURCE_DIR} )
set(CIMG_LIBRARIES ${DEPENDENCIES_LIBRARIES})
