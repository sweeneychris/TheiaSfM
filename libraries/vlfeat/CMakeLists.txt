project(vlfeat)

INCLUDE("${CMAKE_MODULE_PATH}/OptimizeTheiaCompilerFlags.cmake")
OptimizeTheiaCompilerFlags()

INCLUDE_DIRECTORIES(./vl)
SET(vl_sources
  vl/generic.c
  vl/host.c
  vl/imopv.c
  vl/imopv_sse2.c
  vl/mathop.c
  vl/mathop_sse2.c
  vl/random.c
  vl/sift.c)
SET_SOURCE_FILES_PROPERTIES(${vl_sources} PROPERTIES LANGUAGE C)
ADD_LIBRARY(vlfeat SHARED ${vl_sources})
INSTALL(TARGETS vlfeat
  EXPORT  TheiaExport
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib${LIB_SUFFIX}
  ARCHIVE DESTINATION lib${LIB_SUFFIX})
