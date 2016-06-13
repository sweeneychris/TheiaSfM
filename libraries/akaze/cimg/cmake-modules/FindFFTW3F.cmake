# - Try to find FFTW
# Once done this will define
#  FFTW_FOUND - System has FFTW3
#  FFTW_INCLUDE_DIRS - The FFTW3 include directories
#  FFTW_LIBRARIES - The libraries needed to use FFTW3
#  FFTW_DEFINITIONS - Compiler switches required for using FFTW3

find_package(PkgConfig)

pkg_check_modules(PC_FFTW3F QUIET fftw3f)

set(FFTW3F_DEFINITIONS ${PC_FFTW3F_CFLAGS_OTHER})

find_path(FFTW3F_INCLUDE_DIR fftw3.h
          HINTS ${PC_FFTW3F_INCLUDEDIR} ${PC_FFTW3F_INCLUDE_DIRS}
          NAMES fftw3.h )

find_library(FFTW3F_LIBRARY NAMES fftw3f libfftw3f
             HINTS ${PC_FFTW3F_LIBDIR} ${PC_FFTW3F_LIBRARY_DIRS} )

set(FFTW3F_LIBRARIES ${FFTW3F_LIBRARY} )
set(FFTW3F_INCLUDE_DIRS ${FFTW3F_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set FFTW3F_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(FFTW3F  DEFAULT_MSG
                                  FFTW3F_LIBRARY FFTW3F_INCLUDE_DIR)

mark_as_advanced(FFTW3F_INCLUDE_DIR FFTW3F_LIBRARY )