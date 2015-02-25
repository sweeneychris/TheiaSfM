# Copyright (C) 2013 The Regents of the University of California (Regents).
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
#
#     * Neither the name of The Regents or University of California nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#
# Please contact the author of this library if you have any questions.
# Author: Chris Sweeney (cmsweeney@cs.ucsb.edu)

# Add headers
FILE(GLOB THEIA_HDRS *.h)

# Add sources
SET(THEIA_SRC
  image/descriptor/brief_descriptor.cc
  image/descriptor/brisk_descriptor.cc
  image/descriptor/descriptor_extractor.cc
  image/descriptor/freak_descriptor.cc
  image/descriptor/sift_descriptor.cc
  image/image_canvas.cc
  image/keypoint_detector/agast_detector.cc
  image/keypoint_detector/brisk_detector.cc
  image/keypoint_detector/brisk_impl.cc
  image/keypoint_detector/sift_detector.cc
  io/read_bundler_files.cc
  io/read_calibration.cc
  io/read_matches.cc
  io/reconstruction_reader.cc
  io/reconstruction_writer.cc
  io/sift_binary_file.cc
  io/sift_text_file.cc
  io/sift_vsfm_file.cc
  io/write_matches.cc
  math/closed_form_polynomial_solver.cc
  math/polynomial.cc
  math/probability/sequential_probability_ratio.cc
  math/sturm_chain.cc
  util/filesystem.cc
  util/random.cc
  util/stringprintf.cc
  util/threadpool.cc
  util/timer.cc
  matching/cascade_hasher.cc
  matching/cascade_hashing_feature_matcher.cc
  matching/feature_matcher_utils.cc
  sfm/bundle_adjustment/bundle_adjust_two_views.cc
  sfm/bundle_adjustment/bundle_adjustment.cc
  sfm/bundle_adjustment/bundle_adjust_track.cc
  sfm/bundle_adjustment/optimize_relative_position_with_known_rotation.cc
  sfm/camera/camera.cc
  sfm/camera/projection_matrix_utils.cc
  sfm/camera/radial_distortion.cc
  sfm/estimate_track.cc
  sfm/estimate_twoview_info.cc
  sfm/estimators/essential_matrix_estimator.cc
  sfm/estimators/estimate_relative_pose.cc
  sfm/estimators/estimate_triangulation.cc
  sfm/estimators/estimate_uncalibrated_relative_pose.cc
  sfm/estimators/fundamental_matrix_estimator.cc
  sfm/estimators/relative_pose_estimator.cc
  sfm/estimators/uncalibrated_relative_pose_estimator.cc
  sfm/exif_reader.cc
  sfm/feature_extractor.cc
  sfm/filter_view_pairs_from_cycles.cc
  sfm/filter_view_pairs_from_orientation.cc
  sfm/filter_view_pairs_from_relative_translation.cc
  sfm/match_and_verify_features.cc
  sfm/nonlinear_reconstruction_estimator.cc
  sfm/pose/dls_impl.cc
  sfm/pose/dls_pnp.cc
  sfm/pose/eight_point_fundamental_matrix.cc
  sfm/pose/essential_matrix_utils.cc
  sfm/pose/estimate_positions_nonlinear.cc
  sfm/pose/estimate_rotations_nonlinear.cc
  sfm/pose/five_point_focal_length_radial_distortion.cc
  sfm/pose/five_point_relative_pose.cc
  sfm/pose/four_point_focal_length.cc
  sfm/pose/four_point_focal_length_helper.cc
  sfm/pose/four_point_homography.cc
  sfm/pose/four_point_relative_pose_partial_rotation.cc
  sfm/pose/fundamental_matrix_util.cc
  sfm/pose/pairwise_rotation_error.cc
  sfm/pose/pairwise_translation_error.cc
  sfm/pose/perspective_three_point.cc
  sfm/pose/test_util.cc
  sfm/pose/three_point_relative_pose_partial_rotation.cc
  sfm/pose/util.cc
  sfm/reconstruction.cc
  sfm/reconstruction_builder.cc
  sfm/reconstruction_estimator.cc
  sfm/reconstruction_estimator_utils.cc
  sfm/track.cc
  sfm/track_builder.cc
  sfm/transformation/align_point_clouds.cc
  sfm/transformation/gdls_similarity_transform.cc
  sfm/triangulation/triangulation.cc
  sfm/twoview_info.cc
  sfm/verify_two_view_matches.cc
  sfm/view.cc
  sfm/view_graph/orientations_from_view_graph.cc
  sfm/view_graph/remove_disconnected_view_pairs.cc
  sfm/view_graph/triplet_extractor.cc
  sfm/view_graph/view_graph.cc
  )

SET(THEIA_LIBRARY_DEPENDENCIES
  ${CERES_LIBRARIES}
  ${GFLAGS_LIBRARIES}
  ${GLOG_LIBRARIES}
  ${CIMG_LIBRARIES}
  easyexif
  statx
  stlplus3
  vlfeat
  agast)

SET(THEIA_LIBRARY_SOURCE
  ${THEIA_SRC}
  ${THEIA_HDRS})

ADD_LIBRARY(theia ${THEIA_LIBRARY_SOURCE})
SET_TARGET_PROPERTIES(theia PROPERTIES
  VERSION ${THEIA_VERSION}
  SOVERSION ${THEIA_VERSION_MAJOR}
  )
TARGET_LINK_LIBRARIES(theia ${THEIA_LIBRARY_DEPENDENCIES})

INSTALL(TARGETS theia
  EXPORT  TheiaExport
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib${LIB_SUFFIX}
  ARCHIVE DESTINATION lib${LIB_SUFFIX})

IF (BUILD_TESTING)
  MACRO (GTEST FILENAME)
    STRING(REPLACE "/" ";" TEST_NAME_LIST ${FILENAME})
    LIST(REVERSE TEST_NAME_LIST)
    LIST(GET TEST_NAME_LIST 0 TEST_NAME)
    ADD_EXECUTABLE(${TEST_NAME}_test test/test_main.cc ${FILENAME}_test.cc)
    TARGET_LINK_LIBRARIES(${TEST_NAME}_test
      gtest
      theia
      ${THEIA_LIBRARY_DEPENDENCIES})
    ADD_TEST(NAME ${TEST_NAME}_test
      COMMAND ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${TEST_NAME}_test)
  ENDMACRO (GTEST)

  GTEST(image/descriptor/brief_descriptor)
  GTEST(image/descriptor/brisk_descriptor)
  GTEST(image/descriptor/freak_descriptor)
  GTEST(image/descriptor/sift_descriptor)
  GTEST(image/image)
  GTEST(image/keypoint_detector/agast_detector)
  GTEST(image/keypoint_detector/brisk_detector)
  GTEST(image/keypoint_detector/sift_detector)
  GTEST(matching/brute_force_feature_matcher)
  GTEST(matching/distance)
  GTEST(matching/feature_matcher_utils)
  GTEST(math/closed_form_polynomial_solver)
  GTEST(math/graph/connected_components)
  GTEST(math/l1_solver)
  GTEST(math/matrix/gauss_jordan)
  GTEST(math/matrix/rq_decomposition)
  GTEST(math/polynomial)
  GTEST(math/probability/sprt)
  GTEST(sfm/bundle_adjustment/optimize_relative_position_with_known_rotation)
  GTEST(sfm/camera/camera)
  GTEST(sfm/camera/projection_matrix_utils)
  GTEST(sfm/camera/radial_distortion)
  GTEST(sfm/estimators/essential_matrix_estimator)
  GTEST(sfm/estimators/estimate_relative_pose)
  GTEST(sfm/estimators/estimate_triangulation)
  GTEST(sfm/estimators/estimate_uncalibrated_relative_pose)
  GTEST(sfm/estimators/fundamental_matrix_estimator)
  GTEST(sfm/estimators/relative_pose_estimator)
  GTEST(sfm/estimators/uncalibrated_relative_pose_estimator)
  GTEST(sfm/exif_reader)
  GTEST(sfm/filter_view_pairs_from_cycles)
  GTEST(sfm/filter_view_pairs_from_orientation)
  GTEST(sfm/filter_view_pairs_from_relative_translation)
  GTEST(sfm/pose/dls_pnp)
  GTEST(sfm/pose/eight_point_fundamental_matrix)
  GTEST(sfm/pose/essential_matrix_utils)
  GTEST(sfm/pose/estimate_positions_nonlinear)
  GTEST(sfm/pose/five_point_focal_length_radial_distortion)
  GTEST(sfm/pose/five_point_relative_pose)
  GTEST(sfm/pose/four_point_focal_length)
  GTEST(sfm/pose/four_point_homography)
  GTEST(sfm/pose/four_point_relative_pose_partial_rotation)
  GTEST(sfm/pose/fundamental_matrix_util)
  GTEST(sfm/pose/pairwise_rotation_error)
  GTEST(sfm/pose/pairwise_translation_error)
  GTEST(sfm/pose/perspective_three_point)
  GTEST(sfm/pose/three_point_relative_pose_partial_rotation)
  GTEST(sfm/reconstruction)
  GTEST(sfm/track)
  GTEST(sfm/track_builder)
  GTEST(sfm/transformation/align_point_clouds)
  GTEST(sfm/transformation/gdls_similarity_transform)
  GTEST(sfm/triangulation/triangulation)
  GTEST(sfm/twoview_info)
  GTEST(sfm/view)
  GTEST(sfm/view_graph/orientations_from_view_graph)
  GTEST(sfm/view_graph/remove_disconnected_view_pairs)
  GTEST(sfm/view_graph/triplet_extractor)
  GTEST(sfm/view_graph/view_graph)
  GTEST(solvers/arrsac)
  GTEST(solvers/evsac)
  GTEST(solvers/prosac)
  GTEST(solvers/random_sampler)
  GTEST(solvers/ransac)

ENDIF (BUILD_TESTING)
