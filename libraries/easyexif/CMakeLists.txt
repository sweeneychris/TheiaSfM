project(easyexif)

INCLUDE_DIRECTORIES(./)
SET(easyexif_sources
  exif.cpp)

ADD_LIBRARY(easyexif SHARED ${easyexif_sources})
INSTALL(TARGETS easyexif
  EXPORT  TheiaExport
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib${LIB_SUFFIX}
  ARCHIVE DESTINATION lib${LIB_SUFFIX})
