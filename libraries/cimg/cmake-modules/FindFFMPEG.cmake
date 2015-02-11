# - Try to find FFMPEG
# Once done this will define
# 
#  FFMPEG_FOUND          - system has FFMPEG
#  FFMPEG_INCLUDE_DIR    - the include directories
#  FFMPEG_LIBRARY_DIR    - the directory containing the libraries
#  FFMPEG_LIBRARIES      - link these to use FFMPEG
#  FFMPEG_SWSCALE_FOUND  - FFMPEG also has SWSCALE
#   

SET( FFMPEG_HEADERS avformat.h avcodec.h avutil.h )
SET( FFMPEG_PATH_SUFFIXES libavformat libavcodec libavutil )
SET( FFMPEG_SWS_HEADERS swscale.h )
SET( FFMPEG_SWS_PATH_SUFFIXES libswscale )

SET( FFMPEG_LIBRARIES avformat avcodec avutil )
SET( FFMPEG_SWS_LIBRARIES swscale )
INCLUDE(FindPkgConfig)
if ( PKG_CONFIG_FOUND )
    pkg_check_modules( AVFORMAT libavformat )
    pkg_check_modules( AVCODEC libavcodec )
    pkg_check_modules( AVUTIL libavutil )
    pkg_check_modules( SWSCALE libswscale )
endif ( PKG_CONFIG_FOUND )

SET( FFMPEG_LIBRARY_DIR   ${AVFORMAT_LIBRARY_DIRS}
                            ${AVCODEC_LIBRARY_DIRS}
                            ${AVUTIL_LIBRARY_DIRS} )
SET( FFMPEG_INCLUDE_PATHS ${AVFORMAT_INCLUDE_DIRS}
                            ${AVCODEC_INCLUDE_DIRS}
                            ${AVUTIL_INCLUDE_DIRS} )

# add in swscale if found
IF ( SWSCALE_FOUND )
   SET( FFMPEG_LIBRARY_DIR   ${FFMPEG_LIBRARY_DIR}
                             ${SWSCALE_LIBRARY_DIRS} )
   SET( FFMPEG_INCLUDE_PATHS ${FFMPEG_INCLUDE_PATHS}
                             ${SWSCALE_INCLUDE_DIRS} )
   SET( FFMPEG_HEADERS       ${FFMPEG_HEADERS}
                             ${FFMPEG_SWS_HEADERS} )
   SET( FFMPEG_PATH_SUFFIXES ${FFMPEG_PATH_SUFFIXES}
                             ${FFMPEG_SWS_PATH_SUFFIXES} )
   SET( FFMPEG_LIBRARIES     ${FFMPEG_LIBRARIES}
                             ${FFMPEG_SWS_LIBRARIES} )
ENDIF ( SWSCALE_FOUND )

# find includes
SET( INC_SUCCESS 0 )
SET( TMP_ TMP-NOTFOUND )
SET( FFMPEG_INCLUDE_DIR ${FFMPEG_INCLUDE_PATHS} )

FOREACH( INC_ ${FFMPEG_HEADERS} )
#   message( "checking: " ${INC_} )

   FIND_PATH( TMP_ ${INC_}
              PATHS ${FFMPEG_INCLUDE_PATHS}
              PATH_SUFFIXES ${FFMPEG_PATH_SUFFIXES} )
   IF ( TMP_ )
#      message( "found: " ${TMP_} )
      MATH( EXPR INC_SUCCESS ${INC_SUCCESS}+1 )
      SET( FFMPEG_INCLUDE_DIR ${FFMPEG_INCLUDE_DIR} ${TMP_} )
   ENDIF ( TMP_ )
   SET( TMP_ TMP-NOTFOUND )
ENDFOREACH( INC_ )

# clear out duplicates
#LIST( REMOVE_DUPLICATES FFMPEG_INCLUDE_DIR )
#LIST( REMOVE_DUPLICATES FFMPEG_LIBRARY_DIR )

# find the full paths of the libraries
SET( TMP_ TMP-NOTFOUND )
FOREACH( LIB_ ${FFMPEG_LIBRARIES} )
    FIND_LIBRARY( TMP_ NAMES ${LIB_} PATHS ${FFMPEG_LIBRARY_DIR} )
    IF ( TMP_ )
        SET( FFMPEG_LIBRARIES_FULL ${FFMPEG_LIBRARIES_FULL} ${TMP_} )
    ENDIF ( TMP_ )
    SET( TMP_ TMP-NOTFOUND )
ENDFOREACH( LIB_ )
SET ( FFMPEG_LIBRARIES ${FFMPEG_LIBRARIES_FULL} )

LIST( LENGTH FFMPEG_HEADERS LIST_SIZE_ )

SET( FFMPEG_FOUND FALSE )
SET( FFMPEG_SWSCALE_FOUND FALSE )
IF ( ${INC_SUCCESS} EQUAL ${LIST_SIZE_} )
   SET( FFMPEG_FOUND TRUE )
   SET( FFMPEG_SWSCALE_FOUND ${SWSCALE_FOUND} )
ENDIF ( ${INC_SUCCESS} EQUAL ${LIST_SIZE_} ) 
