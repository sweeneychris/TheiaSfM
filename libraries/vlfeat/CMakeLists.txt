project(vlfeat)

include("${CMAKE_MODULE_PATH}/OptimizeTheiaCompilerFlags.cmake")
optimizetheiacompilerflags()

include_directories(./vl)
set(vl_sources
  vl/generic.c
  vl/host.c
  vl/imopv.c
  vl/imopv_sse2.c
  vl/mathop.c
  vl/mathop_sse2.c
  vl/random.c
  vl/sift.c)
set_source_files_properties(${vl_sources} PROPERTIES LANGUAGE C)

if (MSVC)
  add_definitions(-DVL_BUILD_DLL)
  add_definitions(-DVL_DISABLE_SSE2)
endif (MSVC)

add_library(vlfeat SHARED ${vl_sources})
install(TARGETS vlfeat
  EXPORT  TheiaExport
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib${LIB_SUFFIX}
  ARCHIVE DESTINATION lib${LIB_SUFFIX})
