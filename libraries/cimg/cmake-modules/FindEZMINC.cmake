# FindEZMINC.cmake module

# Find the native MINC includes and library
# This module defines
#  EZMINC_INCLUDE_DIR, where to find jpeglib.h, etc.
#  EZMINC_LIBRARIES, the libraries needed to use EZMINC.
#  EZMINC_FOUND, If false, do not try to use EZMINC.
# also defined, but not for general use are
#  EZMINC_LIBRARY, where to find the MINC library.


FIND_PATH(EZMINC_INCLUDE_DIR minc_1_rw.h /usr/include /usr/local/include /usr/local/bic/include /opt/minc-toolkit/include)
FIND_LIBRARY(EZMINC_minc_io_LIBRARY NAMES minc_io HINTS /usr/lib /usr/local/lib /usr/local/bic/lib /opt/minc-toolkit/lib)


IF (EZMINC_INCLUDE_DIR AND EZMINC_minc_io_LIBRARY)
   set(EZMINC_LIBRARIES 
     ${EZMINC_minc_io_LIBRARY} 
     )
   SET(EZMINC_FOUND TRUE)
   
ENDIF (EZMINC_INCLUDE_DIR AND EZMINC_minc_io_LIBRARY)


IF (EZMINC_FOUND)
   IF (NOT Minc_FIND_QUIETLY)
      MESSAGE(STATUS "Found EZMINC: ${EZMINC_LIBRARIES}")
   ENDIF (NOT Minc_FIND_QUIETLY)
ELSE (EZMINC_FOUND)
   IF (Minc_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR "Cound not find EZMINC")
   ENDIF (Minc_FIND_REQUIRED)
ENDIF (EZMINC_FOUND)


mark_as_advanced(
  EZMINC_minc_io_LIBRARY
)
