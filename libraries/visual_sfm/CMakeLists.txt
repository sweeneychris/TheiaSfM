project(visual_sfm)

include_directories(./)
set(visual_sfm_sources
  FeaturePoints.cpp
  MatchFile.cpp)

add_library(visual_sfm STATIC ${visual_sfm_sources})
if (NOT MSVC)
  add_definitions("-Wno-unused-result")
endif (NOT MSVC)
install(TARGETS visual_sfm
  EXPORT  TheiaExport
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib${LIB_SUFFIX}
  ARCHIVE DESTINATION lib${LIB_SUFFIX})
