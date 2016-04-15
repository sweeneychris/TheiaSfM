# Basic programs
add_executable(match_features match_features.cc)
target_link_libraries(match_features theia ${GFLAGS_LIBRARIES})

add_executable(extract_features extract_features.cc)
target_link_libraries(extract_features theia ${GFLAGS_LIBRARIES})

add_executable(compute_two_view_geometry compute_two_view_geometry.cc)
target_link_libraries(compute_two_view_geometry theia ${GFLAGS_LIBRARIES} ${GLOG_LIBRARIES})

# Tools for building SfM models
add_executable(build_reconstruction build_reconstruction.cc)
target_link_libraries(build_reconstruction theia ${GFLAGS_LIBRARIES} ${GLOG_LIBRARIES})

add_executable(build_1dsfm_reconstruction build_1dsfm_reconstruction.cc)
target_link_libraries(build_1dsfm_reconstruction theia ${GFLAGS_LIBRARIES} ${GLOG_LIBRARIES})

# Reconstruction visualization.
find_package(OpenGL REQUIRED)
if (NOT MSVC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-declarations")
  find_package(GLUT REQUIRED)
  include_directories( ${OPENGL_INCLUDE_DIRS}  ${GLUT_INCLUDE_DIRS} )
  set(EXTRA_GL_LIBRARIES ${GLUT_LIBRARIES})
else (NOT MSVC)
  find_package(GLEW REQUIRED)
  include_directories( ${OPENGL_INCLUDE_DIRS}  ${GLEW_INCLUDE_DIRS} )
  set(EXTRA_GL_LIBRARIES ${GLEW_LIBRARIES})
endif (NOT MSVC)

add_executable(view_reconstruction view_reconstruction.cc)
target_link_libraries(view_reconstruction theia ${GFLAGS_LIBRARIES}
  ${OPENGL_LIBRARIES}
  ${EXTRA_GL_LIBRARIES})

# Useful tools for analyzing reconstructions.
add_executable(compute_reconstruction_statistics compute_reconstruction_statistics.cc)
target_link_libraries(compute_reconstruction_statistics theia ${GFLAGS_LIBRARIES} ${GLOG_LIBRARIES})

add_executable(colorize_reconstruction colorize_reconstruction.cc)
target_link_libraries(colorize_reconstruction theia ${GFLAGS_LIBRARIES} ${GLOG_LIBRARIES})

add_executable(compare_reconstructions compare_reconstructions.cc)
target_link_libraries(compare_reconstructions theia ${GFLAGS_LIBRARIES} ${GLOG_LIBRARIES})

add_executable(compute_matching_relative_pose_errors compute_matching_relative_pose_errors.cc)
target_link_libraries(compute_matching_relative_pose_errors theia ${GFLAGS_LIBRARIES} ${GLOG_LIBRARIES})

add_executable(evaluate_relative_translation_optimization evaluate_relative_translation_optimization.cc)
target_link_libraries(evaluate_relative_translation_optimization theia ${GFLAGS_LIBRARIES} ${GLOG_LIBRARIES})

add_executable(verify_1dsfm_input verify_1dsfm_input.cc)
target_link_libraries(verify_1dsfm_input theia ${GFLAGS_LIBRARIES} ${GLOG_LIBRARIES})

# File conversions and exporters.
add_executable(convert_sift_key_file convert_sift_key_file.cc)
target_link_libraries(convert_sift_key_file theia ${GFLAGS_LIBRARIES} ${GLOG_LIBRARIES})

add_executable(convert_bundle_file convert_bundle_file.cc)
target_link_libraries(convert_bundle_file theia ${GFLAGS_LIBRARIES} ${GLOG_LIBRARIES})

add_executable(export_to_nvm_file export_to_nvm_file.cc)
target_link_libraries(export_to_nvm_file theia ${GFLAGS_LIBRARIES} ${GLOG_LIBRARIES})

add_executable(convert_nvm_file convert_nvm_file.cc)
target_link_libraries(convert_nvm_file theia ${GFLAGS_LIBRARIES} ${GLOG_LIBRARIES})

add_executable(convert_theia_reconstruction_to_bundler_file convert_theia_reconstruction_to_bundler_file.cc)
target_link_libraries(convert_theia_reconstruction_to_bundler_file theia ${GFLAGS_LIBRARIES} ${GLOG_LIBRARIES})

add_executable(write_reconstruction_ply_file write_reconstruction_ply_file.cc)
target_link_libraries(write_reconstruction_ply_file theia ${GFLAGS_LIBRARIES} ${GLOG_LIBRARIES})

add_executable(export_reconstruction_to_pmvs export_reconstruction_to_pmvs.cc)
target_link_libraries(export_reconstruction_to_pmvs theia ${GFLAGS_LIBRARIES} ${GLOG_LIBRARIES})

add_executable(create_calibration_file_from_exif create_calibration_file_from_exif.cc)
target_link_libraries(create_calibration_file_from_exif theia ${GFLAGS_LIBRARIES} ${GLOG_LIBRARIES})

add_executable(create_reconstruction_from_strecha_dataset create_reconstruction_from_strecha_dataset.cc)
target_link_libraries(create_reconstruction_from_strecha_dataset theia ${GFLAGS_LIBRARIES} ${GLOG_LIBRARIES})

add_executable(upgrade_matches_file upgrade_matches_file.cc)
target_link_libraries(upgrade_matches_file theia ${GFLAGS_LIBRARIES} ${GLOG_LIBRARIES})

add_executable(upgrade_reconstruction_file upgrade_reconstruction_file.cc)
target_link_libraries(upgrade_reconstruction_file theia ${GFLAGS_LIBRARIES} ${GLOG_LIBRARIES})