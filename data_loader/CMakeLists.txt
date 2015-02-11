ADD_EXECUTABLE(convert_sift_key_file convert_sift_key_file.cc)
TARGET_LINK_LIBRARIES(convert_sift_key_file theia ${GFLAGS_LIBRARIES} ${GLOG_LIBRARIES})

ADD_EXECUTABLE(convert_bundle_file convert_bundle_file.cc)
TARGET_LINK_LIBRARIES(convert_bundle_file theia ${GFLAGS_LIBRARIES} ${GLOG_LIBRARIES})

# A demo program that will load the bigsfm dataset and render it in OpenGL.
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
include_directories( ${OPENGL_INCLUDE_DIRS}  ${GLUT_INCLUDE_DIRS} )

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-declarations")
ADD_EXECUTABLE(view_model view_model.cc)
TARGET_LINK_LIBRARIES(view_model theia ${GFLAGS_LIBRARIES}
  ${OPENGL_LIBRARIES}
  ${GLUT_LIBRARY})
