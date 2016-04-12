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
file(GLOB THEIA_HDRS *.h)

# Add sources
set(THEIA_SRC
  image/descriptor/create_descriptor_extractor.cc
  image/descriptor/descriptor_extractor.cc
  image/descriptor/sift_descriptor.cc
  image/image_canvas.cc
  image/keypoint_detector/sift_detector.cc
  io/import_nvm_file.cc
  io/read_1dsfm.cc
  io/read_bundler_files.cc
  io/read_calibration.cc
  io/read_keypoints_and_descriptors.cc
  io/read_matches.cc
  io/read_matches_deprecated.cc
  io/reconstruction_reader.cc
#  io/reconstruction_reader_deprecated.cc
  io/reconstruction_writer.cc
#  io/reconstruction_writer_deprecated.cc
  io/sift_binary_file.cc
  io/sift_text_file.cc
  io/write_bundler_files.cc
  io/write_keypoints_and_descriptors.cc
  io/write_matches.cc
  io/write_matches_deprecated.cc
  io/write_nvm_file.cc
  io/write_ply_file.cc
  matching/cascade_hasher.cc
  matching/cascade_hashing_feature_matcher.cc
  matching/create_feature_matcher.cc
  matching/feature_matcher_utils.cc
  math/closed_form_polynomial_solver.cc
  math/find_polynomial_roots_companion_matrix.cc
  math/find_polynomial_roots_jenkins_traub.cc
  math/matrix/sparse_cholesky_llt.cc
  math/polynomial.cc
  math/probability/sequential_probability_ratio.cc
  math/qp_solver.cc
  sfm/bundle_adjustment/bundle_adjust_two_views.cc
  sfm/bundle_adjustment/bundle_adjustment.cc
  sfm/bundle_adjustment/create_loss_function.cc
  sfm/bundle_adjustment/optimize_relative_position_with_known_rotation.cc
  sfm/camera/camera.cc
  sfm/camera/projection_matrix_utils.cc
  sfm/camera/radial_distortion.cc
  sfm/camera/shared_extrinsics.cc
  sfm/colorize_reconstruction.cc
  sfm/estimate_track.cc
  sfm/estimate_twoview_info.cc
  sfm/estimators/estimate_calibrated_absolute_pose.cc
  sfm/estimators/estimate_dominant_plane_from_points.cc
  sfm/estimators/estimate_essential_matrix.cc
  sfm/estimators/estimate_fundamental_matrix.cc
  sfm/estimators/estimate_homography.cc
  sfm/estimators/estimate_relative_pose.cc
  sfm/estimators/estimate_similarity_transformation_2d_3d.cc
  sfm/estimators/estimate_triangulation.cc
  sfm/estimators/estimate_uncalibrated_absolute_pose.cc
  sfm/estimators/estimate_uncalibrated_relative_pose.cc
  sfm/exif_reader.cc
  sfm/extract_maximally_parallel_rigid_subgraph.cc
  sfm/feature_extractor.cc
  sfm/feature_extractor_and_matcher.cc
  sfm/filter_view_graph_cycles_by_rotation.cc
  sfm/filter_view_pairs_from_orientation.cc
  sfm/filter_view_pairs_from_relative_translation.cc
  sfm/find_common_tracks_in_views.cc
  sfm/find_common_views_by_name.cc
  sfm/global_pose_estimation/compute_triplet_baseline_ratios.cc
  sfm/global_pose_estimation/least_unsquared_deviation_position_estimator.cc
  sfm/global_pose_estimation/linear_position_estimator.cc
  sfm/global_pose_estimation/linear_rotation_estimator.cc
  sfm/global_pose_estimation/nonlinear_position_estimator.cc
  sfm/global_pose_estimation/nonlinear_rotation_estimator.cc
  sfm/global_pose_estimation/pairwise_rotation_error.cc
  sfm/global_pose_estimation/pairwise_translation_and_scale_error.cc
  sfm/global_pose_estimation/pairwise_translation_error.cc
  sfm/global_pose_estimation/robust_rotation_estimator.cc
  sfm/global_reconstruction_estimator.cc
  sfm/incremental_reconstruction_estimator.cc
  sfm/localize_view_to_reconstruction.cc
  sfm/pose/dls_impl.cc
  sfm/pose/dls_pnp.cc
  sfm/pose/eight_point_fundamental_matrix.cc
  sfm/pose/essential_matrix_utils.cc
  sfm/pose/five_point_focal_length_radial_distortion.cc
  sfm/pose/five_point_relative_pose.cc
  sfm/pose/four_point_focal_length.cc
  sfm/pose/four_point_focal_length_helper.cc
  sfm/pose/four_point_homography.cc
  sfm/pose/four_point_relative_pose_partial_rotation.cc
  sfm/pose/fundamental_matrix_util.cc
  sfm/pose/perspective_three_point.cc
  sfm/pose/seven_point_fundamental_matrix.cc
  sfm/pose/sim_transform_partial_rotation.cc
  sfm/pose/test_util.cc
  sfm/pose/three_point_relative_pose_partial_rotation.cc
  sfm/pose/two_point_pose_partial_rotation.cc
  sfm/pose/util.cc
  sfm/reconstruction.cc
  sfm/reconstruction_builder.cc
  sfm/reconstruction_estimator.cc
  sfm/reconstruction_estimator_utils.cc
  sfm/set_camera_intrinsics_from_priors.cc
  sfm/track.cc
  sfm/track_builder.cc
  sfm/transformation/align_point_clouds.cc
  sfm/transformation/align_reconstructions.cc
  sfm/transformation/align_rotations.cc
  sfm/transformation/gdls_similarity_transform.cc
  sfm/transformation/transform_reconstruction.cc
  sfm/triangulation/triangulation.cc
  sfm/twoview_info.cc
  sfm/verify_two_view_matches.cc
  sfm/view.cc
  sfm/view_graph/orientations_from_maximum_spanning_tree.cc
  sfm/view_graph/remove_disconnected_view_pairs.cc
  sfm/view_graph/triplet_extractor.cc
  sfm/view_graph/view_graph.cc
  util/filesystem.cc
  util/random.cc
  util/stringprintf.cc
  util/threadpool.cc
  util/timer.cc
  )

set(THEIA_LIBRARY_DEPENDENCIES
  ${CERES_LIBRARIES}
  ${GFLAGS_LIBRARIES}
  ${GLOG_LIBRARIES}
  ${CIMG_LIBRARIES}
  ${SUITESPARSE_LIBRARIES}
  easyexif
  statx
  stlplus3
  vlfeat
  visual_sfm)

set(THEIA_LIBRARY_SOURCE
  ${THEIA_SRC}
  ${THEIA_HDRS})

add_library(theia ${THEIA_LIBRARY_SOURCE})
set_target_properties(theia PROPERTIES
  VERSION ${THEIA_VERSION}
  SOVERSION ${THEIA_VERSION_MAJOR}
  )
target_link_libraries(theia ${THEIA_LIBRARY_DEPENDENCIES})

install(TARGETS theia
  EXPORT  TheiaExport
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib${LIB_SUFFIX}
  ARCHIVE DESTINATION lib${LIB_SUFFIX})

if (BUILD_TESTING)
  macro (GTEST FILENAME)
    string(REPLACE "/" ";" TEST_NAME_LIST ${FILENAME})
    list(REVERSE TEST_NAME_LIST)
    list(GET TEST_NAME_LIST 0 TEST_NAME)
    add_executable(${TEST_NAME}_test test/test_main.cc ${FILENAME}_test.cc)
    target_link_libraries(${TEST_NAME}_test
      gtest
      theia
      ${THEIA_LIBRARY_DEPENDENCIES})
    add_test(NAME ${TEST_NAME}_test
      COMMAND ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${TEST_NAME}_test)
  endmacro (GTEST)

  gtest(image/descriptor/sift_descriptor)
  gtest(image/image)
  gtest(image/keypoint_detector/sift_detector)
  gtest(matching/brute_force_feature_matcher)
  gtest(matching/cascade_hashing_feature_matcher)
  gtest(matching/distance)
  gtest(matching/feature_matcher_utils)
  gtest(math/closed_form_polynomial_solver)
  gtest(math/find_polynomial_roots_companion_matrix)
  gtest(math/find_polynomial_roots_jenkins_traub)
  gtest(math/graph/connected_components)
  gtest(math/graph/minimum_spanning_tree)
  gtest(math/graph/normalized_graph_cut)
  gtest(math/l1_solver)
  gtest(math/matrix/gauss_jordan)
  gtest(math/matrix/rq_decomposition)
  gtest(math/polynomial)
  gtest(math/probability/sprt)
  gtest(math/qp_solver)
  gtest(sfm/bundle_adjustment/optimize_relative_position_with_known_rotation)
  gtest(sfm/camera/camera)
  gtest(sfm/camera/projection_matrix_utils)
  gtest(sfm/camera/radial_distortion)
  gtest(sfm/estimators/estimate_calibrated_absolute_pose)
  gtest(sfm/estimators/estimate_dominant_plane_from_points)
  gtest(sfm/estimators/estimate_essential_matrix)
  gtest(sfm/estimators/estimate_fundamental_matrix)
  gtest(sfm/estimators/estimate_homography)
  gtest(sfm/estimators/estimate_relative_pose)
  gtest(sfm/estimators/estimate_similarity_transformation_2d_3d)
  gtest(sfm/estimators/estimate_triangulation)
  gtest(sfm/estimators/estimate_uncalibrated_absolute_pose)
  gtest(sfm/estimators/estimate_uncalibrated_relative_pose)
  gtest(sfm/exif_reader)
  gtest(sfm/extract_maximally_parallel_rigid_subgraph)
  gtest(sfm/filter_view_graph_cycles_by_rotation)
  gtest(sfm/filter_view_pairs_from_orientation)
  gtest(sfm/filter_view_pairs_from_relative_translation)
  gtest(sfm/find_common_tracks_in_views)
  gtest(sfm/find_common_views_by_name)
  gtest(sfm/global_pose_estimation/compute_triplet_baseline_ratios)
  gtest(sfm/global_pose_estimation/least_unsquared_deviation_position_estimator)
  gtest(sfm/global_pose_estimation/linear_position_estimator)
  gtest(sfm/global_pose_estimation/linear_rotation_estimator)
  gtest(sfm/global_pose_estimation/nonlinear_position_estimator)
  gtest(sfm/global_pose_estimation/pairwise_rotation_error)
  gtest(sfm/global_pose_estimation/pairwise_translation_and_scale_error)
  gtest(sfm/global_pose_estimation/pairwise_translation_error)
  gtest(sfm/global_pose_estimation/robust_rotation_estimator)
  gtest(sfm/pose/dls_pnp)
  gtest(sfm/pose/eight_point_fundamental_matrix)
  gtest(sfm/pose/essential_matrix_utils)
  gtest(sfm/pose/five_point_focal_length_radial_distortion)
  gtest(sfm/pose/five_point_relative_pose)
  gtest(sfm/pose/four_point_focal_length)
  gtest(sfm/pose/four_point_homography)
  gtest(sfm/pose/four_point_relative_pose_partial_rotation)
  gtest(sfm/pose/fundamental_matrix_util)
  gtest(sfm/pose/perspective_three_point)
  gtest(sfm/pose/seven_point_fundamental_matrix)
  gtest(sfm/pose/sim_transform_partial_rotation)
  gtest(sfm/pose/three_point_relative_pose_partial_rotation)
  gtest(sfm/pose/two_point_pose_partial_rotation)
  gtest(sfm/reconstruction)
  gtest(sfm/track)
  gtest(sfm/track_builder)
  gtest(sfm/transformation/align_point_clouds)
  gtest(sfm/transformation/align_reconstructions)
  gtest(sfm/transformation/align_rotations)
  gtest(sfm/transformation/gdls_similarity_transform)
  gtest(sfm/triangulation/triangulation)
  gtest(sfm/twoview_info)
  gtest(sfm/view)
  gtest(sfm/view_graph/orientations_from_maximum_spanning_tree)
  gtest(sfm/view_graph/remove_disconnected_view_pairs)
  gtest(sfm/view_graph/triplet_extractor)
  gtest(sfm/view_graph/view_graph)
  gtest(solvers/arrsac)
  gtest(solvers/evsac)
  gtest(solvers/lmed)
  gtest(solvers/prosac)
  gtest(solvers/random_sampler)
  gtest(solvers/ransac)
  gtest(util/mutable_priority_queue)
  gtest(util/lru_cache)
endif (BUILD_TESTING)
