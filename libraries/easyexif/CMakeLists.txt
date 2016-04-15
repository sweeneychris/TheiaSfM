project(easyexif)

include_directories(./)

if (NOT MSVC)
  include(CheckCXXCompilerFlag)
  check_cxx_compiler_flag("-std=c++11" COMPILER_HAS_CXX11_FLAG)
  if (COMPILER_HAS_CXX11_FLAG)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
    if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
      # Mac OS X before Mavericks uses libstdc++ by default but does not support
      # c++11. Force it to use libc++.
      if (CMAKE_SYSTEM_NAME MATCHES "Darwin")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
      endif (CMAKE_SYSTEM_NAME MATCHES "Darwin")
    endif (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  else (COMPILER_HAS_CXX11_FLAG)
    message(FATAL_ERROR "A compiler with C++11 support is required for EasyExif.")
  endif (COMPILER_HAS_CXX11_FLAG)
endif (NOT MSVC)


set(easyexif_sources
  exif.cpp)

add_library(easyexif STATIC ${easyexif_sources})
install(TARGETS easyexif
  EXPORT  TheiaExport
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib${LIB_SUFFIX}
  ARCHIVE DESTINATION lib${LIB_SUFFIX})
