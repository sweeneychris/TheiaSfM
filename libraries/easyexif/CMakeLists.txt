project(easyexif)

include_directories(./)
set(easyexif_sources
  exif.cpp)

add_library(easyexif SHARED ${easyexif_sources})
install(TARGETS easyexif
  EXPORT  TheiaExport
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib${LIB_SUFFIX}
  ARCHIVE DESTINATION lib${LIB_SUFFIX})
