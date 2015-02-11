project(stlplus3)

INCLUDE_DIRECTORIES(./)
SET(stlplus3_sources
  file_system.cpp
  portability_fixes.cpp
  wildcard.cpp)

ADD_LIBRARY(stlplus3 SHARED ${stlplus3_sources})
INSTALL(TARGETS stlplus3
  EXPORT  TheiaExport
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib${LIB_SUFFIX}
  ARCHIVE DESTINATION lib${LIB_SUFFIX})
