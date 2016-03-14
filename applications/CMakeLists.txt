# Basic programs
ADD_EXECUTABLE(match_features match_features.cc)
TARGET_LINK_LIBRARIES(match_features theia ${GFLAGS_LIBRARIES})

ADD_EXECUTABLE(extract_features extract_features.cc)
TARGET_LINK_LIBRARIES(extract_features theia ${GFLAGS_LIBRARIES})

ADD_EXECUTABLE(compute_two_view_geometry compute_two_view_geometry.cc)
TARGET_LINK_LIBRARIES(compute_two_view_geometry theia ${GFLAGS_LIBRARIES} ${GLOG_LIBRARIES})

# Tools for building SfM models
ADD_EXECUTABLE(build_reconstruction build_reconstruction.cc)
TARGET_LINK_LIBRARIES(build_reconstruction theia ${GFLAGS_LIBRARIES} ${GLOG_LIBRARIES})

ADD_EXECUTABLE(build_1dsfm_reconstruction build_1dsfm_reconstruction.cc)
TARGET_LINK_LIBRARIES(build_1dsfm_reconstruction theia ${GFLAGS_LIBRARIES} ${GLOG_LIBRARIES})

# Reconstruction visualization.
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
include_directories( ${OPENGL_INCLUDE_DIRS}  ${GLUT_INCLUDE_DIRS} )

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-declarations")
ADD_EXECUTABLE(view_reconstruction view_reconstruction.cc)
TARGET_LINK_LIBRARIES(view_reconstruction theia ${GFLAGS_LIBRARIES}
  ${OPENGL_LIBRARIES}
  ${GLUT_LIBRARY})

# Useful tools for analyzing reconstructions.
ADD_EXECUTABLE(compute_reconstruction_statistics compute_reconstruction_statistics.cc)
TARGET_LINK_LIBRARIES(compute_reconstruction_statistics theia ${GFLAGS_LIBRARIES} ${GLOG_LIBRARIES})

ADD_EXECUTABLE(colorize_reconstruction colorize_reconstruction.cc)
TARGET_LINK_LIBRARIES(colorize_reconstruction theia ${GFLAGS_LIBRARIES} ${GLOG_LIBRARIES})

ADD_EXECUTABLE(compare_reconstructions compare_reconstructions.cc)
TARGET_LINK_LIBRARIES(compare_reconstructions theia ${GFLAGS_LIBRARIES} ${GLOG_LIBRARIES})

ADD_EXECUTABLE(compute_matching_relative_pose_errors compute_matching_relative_pose_errors.cc)
TARGET_LINK_LIBRARIES(compute_matching_relative_pose_errors theia ${GFLAGS_LIBRARIES} ${GLOG_LIBRARIES})

ADD_EXECUTABLE(evaluate_relative_translation_optimization evaluate_relative_translation_optimization.cc)
TARGET_LINK_LIBRARIES(evaluate_relative_translation_optimization theia ${GFLAGS_LIBRARIES} ${GLOG_LIBRARIES})

ADD_EXECUTABLE(verify_1dsfm_input verify_1dsfm_input.cc)
TARGET_LINK_LIBRARIES(verify_1dsfm_input theia ${GFLAGS_LIBRARIES} ${GLOG_LIBRARIES})

# File conversions and exporters.
ADD_EXECUTABLE(convert_sift_key_file convert_sift_key_file.cc)
TARGET_LINK_LIBRARIES(convert_sift_key_file theia ${GFLAGS_LIBRARIES} ${GLOG_LIBRARIES})

ADD_EXECUTABLE(convert_bundle_file convert_bundle_file.cc)
TARGET_LINK_LIBRARIES(convert_bundle_file theia ${GFLAGS_LIBRARIES} ${GLOG_LIBRARIES})

ADD_EXECUTABLE(export_to_nvm_file export_to_nvm_file.cc)
TARGET_LINK_LIBRARIES(export_to_nvm_file theia ${GFLAGS_LIBRARIES} ${GLOG_LIBRARIES})

ADD_EXECUTABLE(convert_nvm_file convert_nvm_file.cc)
TARGET_LINK_LIBRARIES(convert_nvm_file theia ${GFLAGS_LIBRARIES} ${GLOG_LIBRARIES})

ADD_EXECUTABLE(convert_theia_reconstruction_to_bundler_file convert_theia_reconstruction_to_bundler_file.cc)
TARGET_LINK_LIBRARIES(convert_theia_reconstruction_to_bundler_file theia ${GFLAGS_LIBRARIES} ${GLOG_LIBRARIES})

ADD_EXECUTABLE(write_reconstruction_ply_file write_reconstruction_ply_file.cc)
TARGET_LINK_LIBRARIES(write_reconstruction_ply_file theia ${GFLAGS_LIBRARIES} ${GLOG_LIBRARIES})

ADD_EXECUTABLE(export_reconstruction_to_pmvs export_reconstruction_to_pmvs.cc)
TARGET_LINK_LIBRARIES(export_reconstruction_to_pmvs theia ${GFLAGS_LIBRARIES} ${GLOG_LIBRARIES})

ADD_EXECUTABLE(create_calibration_file_from_exif create_calibration_file_from_exif.cc)
TARGET_LINK_LIBRARIES(create_calibration_file_from_exif theia ${GFLAGS_LIBRARIES} ${GLOG_LIBRARIES})

ADD_EXECUTABLE(create_reconstruction_from_strecha_dataset create_reconstruction_from_strecha_dataset.cc)
TARGET_LINK_LIBRARIES(create_reconstruction_from_strecha_dataset theia ${GFLAGS_LIBRARIES} ${GLOG_LIBRARIES})

ADD_EXECUTABLE(upgrade_matches_file upgrade_matches_file.cc)
TARGET_LINK_LIBRARIES(upgrade_matches_file theia ${GFLAGS_LIBRARIES} ${GLOG_LIBRARIES})

ADD_EXECUTABLE(upgrade_reconstruction_file upgrade_reconstruction_file.cc)
TARGET_LINK_LIBRARIES(upgrade_reconstruction_file theia ${GFLAGS_LIBRARIES} ${GLOG_LIBRARIES})