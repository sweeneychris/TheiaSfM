.. highlight:: c++

.. default-domain:: cpp

.. _documentation-pose:

=====================
Pose and Resectioning
=====================

Theia contains efficient and robust implementations of the following pose and
resectioning algorithms. We attempted to make each method as general as possible so that users were not tied to Theia data structures to use the methods. The interface for all pose methods uses Eigen types for feature positions, 3D positions, and pose rotations and translations.

.. _section-p3p:

Perspective Three Point (P3P)
=============================

  .. function:: bool PoseFromThreePoints(const Eigen::Vector2d feature_position[3], const Eigen::Vector3d world_point[3], std::vector<Eigen::Matrix3d>* solution_rotations, std::vector<Eigen::Vector3d>* solution_translations)

    Computes camera pose using the three point algorithm and returns all
    possible solutions (up to 4). Follows steps from the paper "A Novel
    Parameterization of the Perspective-Three-Point Problem for a direct
    computation of Absolute Camera position and Orientation" by [Kneip]_\. This
    algorithm has been proven to be up to an order of magnitude faster than
    other methods. The output rotation and translation define world-to-camera
    transformation.

    ``feature_position``: Image points corresponding to model points. These should be
    calibrated image points as opposed to pixel values.

    ``world_point``: 3D location of features.

    ``solution_rotations``: the rotation matrix of the candidate solutions

    ``solution_translation``: the translation of the candidate solutions

    ``returns``: Whether the pose was computed successfully, along with the
    output parameters ``rotation`` and ``translation`` filled with the valid
    poses.

.. _section-five_point_essential_matrix:

Five Point Relative Pose
========================

  .. function:: bool FivePointRelativePose(const Eigen::Vector2d image1_points[5], const Eigen::Vector2d image2_points[5], std::vector<Eigen::Matrix3d>* rotation, std::vector<Eigen::Vector3d>* translation)

    Computes the relative pose between two cameras using 5 corresponding
    points. Algorithm is implemented based on "Recent Developments on Direct
    Relative Orientation" by [Stewenius5pt]_. This algorithm is known to be more
    numerically stable while only slightly slower than the [Nister]_ method. The
    rotation and translation returned are defined such that
    :math:`E=[t]_{\times} * R` and :math:`y^\top * E * x = 0` where :math:`y`
    are points from image2 and :math:`x` are points from image1.

    ``image1_points``: Location of features on the image plane of image 1.

    ``image2_points``: Location of features on the image plane of image 2.

    ``returns``: Output the number of poses computed as well as the relative
    rotation and translation.


.. _section-four_point_homography:

Four Point Algorithm for Homography
===================================

  .. function:: bool FourPointHomography(const std::vector<Eigen::Vector2d>& image_1_points, const std::vector<Eigen::Vector2d>& image_2_points, Eigen::Matrix3d* homography)

    Computes the 2D `homography
    <http://en.wikipedia.org/wiki/Homography_(computer_vision)>`_ mapping points
    in image 1 to image 2 such that: :math:`x' = Hx` where :math:`x` is a point in
    image 1 and :math:`x'` is a point in image 2. The algorithm implemented is
    the DLT algorithm based on algorithm 4.2 in [HartleyZisserman]_.

    ``image_1_points``: Image points from image 1. At least 4 points must be
    passed in.

    ``image_2_points``: Image points from image 2. At least 4 points must be
    passed in.

    ``homography``: The computed 3x3 homography matrix.

.. _section-eight_point:

Eight Point Algorithm for Fundamental Matrix
============================================

  .. function:: bool EightPointFundamentalMatrix(const std::vector<Eigen::Vector2d>& image_1_points, const std::vector<Eigen::Vector2d>& image_2_points, Eigen::Matrix3d* fundamental_matrix)

    Computes the `fundamental matrix
    <http://en.wikipedia.org/wiki/Fundamental_matrix_(computer_vision)>`_ relating
    image points between two images such that :math:`x' F x = 0` for all
    correspondences :math:`x` and :math:`x'` in images 1 and 2 respectively. The
    normalized eight point algorithm is a speedy estimation of the fundamental
    matrix (Alg 11.1 in [HartleyZisserman]_) that minimizes an algebraic error.

    ``image_1_points``: Image points from image 1. At least 8 points must be
    passed in.

    ``image_2_points``: Image points from image 2. At least 8 points must be
    passed in.

    ``fundamental_matrix``: The computed fundamental matrix.

    ``returns:`` true on success, false on failure.

.. _section-dls_pnp:

Perspective N-Point
===================

  .. function:: void DlsPnp(const std::vector<Eigen::Vector2d>& feature_position, const std::vector<Eigen::Vector3d>& world_point, std::vector<Eigen::Quaterniond>* solution_rotation, std::vector<Eigen::Vector3d>* solution_translation)

    Computes the camera pose using the Perspective N-point method from "A Direct
    Least-Squares (DLS) Method for PnP" by [Hesch]_ and Stergios Roumeliotis. This
    method is extremely scalable and highly accurate for the PnP problem. A
    minimum of 4 points are required, but there is no maximum number of points
    allowed as this is a least-squared approach. Theoretically, up to 27 solutions
    may be returned, but in practice only 4 real solutions arise and in almost all
    cases where n >= 6 there is only one solution which places the observed points
    in front of the camera. The returned rotation and translations are
    world-to-camera transformations.

    ``feature_position``: Normalized image rays corresponding to model points. Must
    contain at least 4 points.

    ``points_3d``: 3D location of features. Must correspond to the image_ray of
    the same index. Must contain the same number of points as image_ray, and at
    least 4.

    ``solution_rotation``: the rotation quaternion of the candidate solutions

    ``solution_translation``: the translation of the candidate solutions


.. _section-four_point_focal_length:

Four Point Focal Length
=======================

  .. function:: int FourPointPoseAndFocalLength(const std::vector<Eigen::Vector2d>& feature_positions, const std::vector<Eigen::Vector3d>& world_points, std::vector<Eigen::Matrix<double, 3, 4> >* projection_matrices)

    Computes the camera pose and unknown focal length of an image given four 2D-3D
    correspondences, following the method of [Bujnak]_. This method involves
    computing a grobner basis from a modified constraint of the focal length and
    pose projection.

    ``feature_position``: Normalized image rays corresponding to model points. Must
    contain at least 4 points.

    ``points_3d``: 3D location of features. Must correspond to the image_ray of
    the same index. Must contain the same number of points as image_ray, and at
    least 4.

    ``projection_matrices``: The solution world-to-camera projection matrices,
    inclusive of the unknown focal length. For a focal length f and a camera
    calibration matrix :math:`K=diag(f, f, 1)`, the projection matrices returned
    are of the form :math:`P = K * [R | t]`.


.. _section-five_point_focal_length_radial_distortion:

Five Point Focal Length and Radial Distortion
=============================================

  .. function:: bool FivePointFocalLengthRadialDistortion(const std::vector<Eigen::Vector2d>& feature_positions, const std::vector<Eigen::Vector3d>& world_points, const int num_radial_distortion_params, std::vector<Eigen::Matrix<double, 3, 4> >* projection_matrices, std::vector<std::vector<double> >* radial_distortions)

    Compute the absolute pose, focal length, and radial distortion of a camera
    using five 3D-to-2D correspondences [Kukelova]_. The method solves for the
    projection matrix (up to scale) by using a cross product constraint on the
    standard projection equation. This allows for simple solution to the first two
    rows of the projection matrix, and the third row (which contains the focal
    length and distortion parameters) can then be solved with SVD on the remaining
    constraint equations from the first row of the projection matrix. See the
    paper for more details.

    ``feature_positions``: the 2D location of image features. Exactly five
    features must be passed in.

    ``world_points``: 3D world points corresponding to the features
    observed. Exactly five points must be passed in.

    ``num_radial_distortion_params``: The number of radial distortion paramters to
	solve for. Must be 1, 2, or 3.

    ``projection_matrices``: Camera projection matrices (that encapsulate focal
	length). These solutions are only valid up to scale.

    ``radial_distortions``: Each entry of this vector contains a vector with the
    radial distortion parameters (up to 3, but however many were specified in
    ``num_radial_distortion_params``).

    ``return``: true if successful, false if not.

Three Point Relative Pose with a Partially Known Rotation
=========================================================

  .. function:: void ThreePointRelativePosePartialRotation(const Eigen::Vector3d& rotation_axis, const Eigen::Vector3d image_1_rays[3], const Eigen::Vector3d image_2_rays[3], std::vector<Eigen::Quaterniond>* soln_rotations, std::vector<Eigen::Vector3d>* soln_translations)

    Computes the relative pose between two cameras using 3 correspondences and a
    known vertical direction as a Quadratic Eigenvalue Problem [SweeneyQEP]_. Up
    to 6 solutions are returned such that :math:`x_2 = R * x_1 + t` for rays
    :math:`x_1` in image 1 and rays :math:`x_2` in image 2. The ``axis`` that is
    passed in as a known axis of rotation (when considering rotations as an
    angle axis). This is equivalent to aligning the two cameras to a common
    direction such as the vertical direction, which can be done using IMU data.

Four Point Relative Pose with a Partially Known Rotation
========================================================

  .. function:: void FourPointRelativePosePartialRotation(const Eigen::Vector3d& rotation_axis, const Eigen::Vector3d image_1_origins[4], const Eigen::Vector3d image_1_rays[4], const Eigen::Vector3d image_2_origins[4], const Eigen::Vector3d image_2_rays[4], std::vector<Eigen::Quaterniond>* soln_rotations, std::vector<Eigen::Vector3d>* soln_translations)

    Computes the relative pose between two generalized cameras using 4
    correspondences and a known vertical direction as a Quadratic Eigenvalue
    Problem [SweeneyQEP]_. A generalized camera is a camera setup with multiple
    cameras such that the cameras do not have the same center of projection
    (e.g., a multi-camera rig mounted on a car). Up to 8 solutions are returned
    such that :math:`x_2 = R * x_1 + t` for rays :math:`x_1` in image 1 and rays
    :math:`x_2` in image 2. The axis that is passed in as a known axis of
    rotation (when considering rotations as an angle axis). This is equivalent
    to aligning the two cameras to a common direction such as the vertical
    direction, which can be done using IMU data.


Two Point Absolute Pose with a Partially Known Rotation
=======================================================

  .. function:: int TwoPointPosePartialRotation(const Eigen::Vector3d& axis, const Eigen::Vector3d& model_point_1, const Eigen::Vector3d& model_point_2, const Eigen::Vector3d& image_ray_1, const Eigen::Vector3d& image_ray_2, Eigen::Quaterniond soln_rotations[2], Eigen::Vector3d soln_translations[2])


    Solves for the limited pose of a camera from two 3D points to image ray
    correspondences. The pose is limited in that while it solves for the three
    translation components, it only solves for a single rotation around a passed
    axis.

    This is intended for use with camera phones that have accelerometers, so that
    the 'up' vector is known, meaning the other two rotations are known. The
    effect of the other rotations should be removed before using this function.

    This implementation is intended to form the core of a RANSAC routine, and as
    such has an optimized interface for this use case.

    Computes the limited pose between the 3D model points and the (unit-norm)
    image rays. Places the rotation and translation solutions in soln_rotations
    and soln_translations.
    There are at most 2 solutions, and the number of solutions is returned.

    The rotations and translation are defined such that model points are
    transformed according to  :math:`image_point = Q * model_point + t`

    This function computes the rotation and translation such that the model
    points, after transformation, lie along the corresponding image_rays. The
    axis referred to is the axis of rotation between the camera coordinate system
    and world (3D point) coordinate system. For most users, this axis will be
    (0, 1, 0) i.e., the up direction. This requires that the input image rays
    have been rotated such that the up direction of the camera coordinate system
    is indeed equal to (0, 1, 0).

    When using this algorithm please cite the paper [SweeneyISMAR]_.
