.. highlight:: c++

.. default-domain:: cpp

.. _documentation-sfm:

===========================
Structure from Motion (SfM)
===========================

Theia has a full Structure-from-Motion pipeline that is extremely efficient. Our
overall pipeline consists of several steps. First, we extract features (SIFT is
the default). Then, we perform two-view matching and geometric verification to
obtain relative poses between image pairs and create a :class:`ViewGraph`. Next,
we perform global pose estimation with gloabl SfM. Global SfM is different from
incremental SfM in that it considers the entire view graph at the same time
instead of incrementally adding more and more images to the
:class:`Reconstruction`. Global SfM methods have been proven to be very fast
with comparable or better accuracy to incremental SfM approaches (See
[JiangICCV]_, [MoulonICCV]_, [WilsonECCV]_), and they are much more readily
parallelized. After we have obtained camera poses, we perform triangulation and
:class:`BundleAdjustment` to obtain a valid 3D reconstruction consisting of
cameras and 3D points.

The first step towards creating a reconstruction is to determine images which
view the same objects. To do this, we must create a :class:`ViewGraph`.

  #. Extract features in images.
  #. Match features to obtain image correspondences.
  #. Estimate camera poses from two-view matches and geometries.

INSERT FIGURE HERE!!!

#1. and #2. have been covered in other sections, so we will focus on creating a
reconstruction from two-view matches and geometry. First, we will describe the
fundamental elements of our reconstruction.

Reconstruction
==============

.. class:: Reconstruction

INSERT FIGURE HERE

At the core of our SfM pipeline is an SfM :class:`Reconstruction`. A
:class:`Reconstruction` is the representation of a 3D reconstuction consisting
of a set of unique :class:`Views` and :class:`Tracks`. More importantly, the
:class:`Reconstruction` class contains visibility information relating all of
the Views and Tracks to each other. We identify each :class:`View` uniquely
based on the name (a string). A good name for the view is the filename of the
image that corresponds to the :class:`View`

When creating an SfM reconstruction, you should add each :class:`View` and
:class:`Track` through the :class:`Reconstruction` object. This will ensure that
visibility information (such as which Tracks are observed a given View and which
Views see a given Track) stays accurate. Views and Tracks are given a unique ID
when added to the :class:`Reconstruction` to help make use of these structures
lightweight and efficient.

.. function:: ViewId AddView(const std::string& view_name)

    Adds a view to the reconstruction with the default initialization. The ViewId
    returned is guaranteed to be unique or will be kInvalidViewId if the method
    fails. Each view is uniquely identified by the view name (a good view name could
    be the filename of the image).

.. function:: bool RemoveView(const ViewId view_id)

    Removes the view from the reconstruction and removes all references to the view in
    the tracks.

    .. NOTE:: Any tracks that have length 0 after the view is removed will also be removed.

.. function:: int NumViews() const
.. function:: int NumTracks() const

.. function:: const View* View(const ViewId view_id) const
.. function:: View* MutableView(const ViewId view_id)

    Returns the View or a nullptr if the track does not exist.

.. function:: std::vector<ViewId> ViewIds() const

    Return all ViewIds in the reconstruction.

.. function:: ViewId ViewIdFromName(const std::string& view_name) const

    Returns to ViewId of the view name, or kInvalidViewId if the view does not
    exist.

.. function:: TrackId AddTrack(const std::vector<std::pair<ViewId, Feature> >& track)

    Add a track to the reconstruction with all of its features across views that observe
    this track. Each pair contains a feature and the corresponding View name
    (i.e., the string) that observes the feature. A new View will be created if
    a View with the view name does not already exist. This method will not
    estimate the position of the track. The TrackId returned will be unique or
    will be kInvalidTrackId if the method fails.

.. function:: bool RemoveTrack(const TrackId track_id)

    Removes the track from the reconstruction and from any Views that observe this
    track. Returns true on success and false on failure (e.g., the track does
    not exist).

.. function:: const Track* Track(const TrackId track_id) const
.. function:: Track* MutableTrack(const TrackId track_id)

    Returns the Track or a nullptr if the track does not exist.

.. function:: std::vector<TrackId> TrackIds() const

    Return all TrackIds in the reconstruction.

ViewGraph
=========

.. class:: ViewGraph

INSERT FIGURE HERE

A :class:`ViewGraph` is a basic SfM construct that is created from two-view
matching information. Any pair of views that have a view correlation form an
edge in the :class:`ViewGraph` such that the nodes in the graph are
:class:`View` that are connected by :class:`TwoViewInfo` objects that contain
information about the relative pose between the Views as well as matching
information.

Once you have a set of views and match information, you can add them to the view graph:

.. code:: c++

  std::vector<View> views;
  // Match all views in the set.
  std::vector<ViewIdPair, TwoViewInfo> view_pair_matches;

  ViewGraph view_graph;
  for (const auto& view_pair : view_pair_matches) {
    const ViewIdPair& view_id_pair = view_pair.first;
    const TwoViewInfo& two_view_info = view_pair.second;
    // Only add view pairs to the view graph if they have strong visual coherence.
    if (two_view_info.num_matched_features > min_num_matched_features) {
      view_graph.AddEdge(views[view_id_pair.first],
                         views[view_id_pair.second],
                         two_view_info);
    }
  }

  // Process and/or manipulate the view graph.

The edge values are especially useful for one-shot SfM where the relative poses
are heavily exploited for computing the final poses. Without a proper
:class:`ViewGraph`, one-shot SfM would not be possible.

Views and Tracks
================

.. class:: View

At the heart of our SfM framework is the :class:`View` class which represents
everything about an image that we want to reconstruct. It contains information
about features from the image, camera pose information, and EXIF
information. Views make up our basic visiblity constraints and are a fundamental
part of the SfM pipeline.

.. class:: Track

A :class:`Track` represents a feature that has been matached over potentially
many images. When a feature appears in multiple images it typically means that
the features correspond to the same 3D point. These 3D points are useful
constraints in SfM reconstruction, as they represent the "structure" in
"Structure-from-Motion" and help to build a point cloud for our reconstruction.

TwoViewInfo
===========

.. class:: TwoViewInfo

After image matching is performed we can create a :class:`ViewGraph` that
explains the relative pose information between two images that have been
matched. The :class:`TwoViewInfo` struct is specified as:

.. code:: c++

  struct TwoViewInfo {
    double focal_length_1;
    double focal_length_2;

    Eigen::Vector3d position_2;
    Eigen::Vector3d rotation_2;

    // Number of features that were matched and geometrically verified betwen the
    // images.
    int num_verified_matches;
  };

This information serves the purpose of an edge in the view graph that describes
visibility information between all views. The relative poses here are used to
estimate global poses for the cameras.

Camera
======

.. class:: Camera

Each :class:`View` contains a :class:`Camera` object that contains intrinsic and
extrinsic information about the camera that observed the scene. Theia has an
efficient, compact :class:`Camera` class that abstracts away common image
operations. This greatly relieves the pain of manually dealing with calibration
and geometric transformations of images. We represent camera intrinsics such
that the calibration matrix is:

.. math::

  K = \left[\begin{matrix}f & s & p_x \\ 0 & f * a & p_y \\ 0 & 0 & 1 \end{matrix} \right]

where :math:`f` is the focal length (in pixels), :math:`s` is the skew,
:math:`a` is the aspect ratio and :math:`p` is the principle point of the
camera. All of these intrinsics may be accessed with getter and setter methods,
e.g., :code:`double GetFocalLenth()` or :code:`void SetFocalLength(const double
focal_length)`. Note that we do additionally allow for up to two radial
distortion parameters, but these are not part of the calibration matrix so they
must be set or retrieved separately from the corresponding getter/setter
methods.

We store the camera pose information as the transformation which maps world
coordinates into camera coordinates. Our rotation is stored internally as an
`SO(3)` rotation, which makes optimization with :class:`BundleAdjustment` more
effective since the value is always a valid rotation (unlike e.g., Quaternions
that must be normalized after each optimization step). However, for convenience
we provide an interface to retrieve the rotation as a rotation matrix as
well. Further, we store the camera position as opposed to the translation.

The convenience of this camera class is clear with the common example of 3D
point reprojection.

.. code:: c++

   // Open an image and obtain camera parameters.
   FloatImage image("my_image.jpg");
   double focal_length;
   CHECK(image.FocalLengthPixels(&focal_length));
   const double radial_distortion1 = value obtained elsewhere...
   const double radial_distortion2 = value obtained elsewhere...
   const Eigen::Matrix3d rotation = value obtained elsewhere...
   const Eigen::Vector3d position = value obtained elsewhere...

   // Set up the camera.
   Camera camera;
   camera.SetOrientationFromRotationMatrix(rotation);
   camera.SetPosition(position);
   camera.SetFocalLength(focal_length);
   camera.SetPrincipalPoint(image.Width() / 2.0, image.Height() / 2.0);
   camera.SetRadialDistortion(radial_distortion1, radial_distortion2);

   // Obtain a homogeneous 3D point
   const Eigen::Vector4d homogeneous_point3d = value obtained elsewhere...

   // Reproject the 3D point to a pixel.
   Eigen::Vector2d reprojection_pixel;
   const double depth = camera.ProjectPoint(homogeneous_point3d, &pixel);
   if (depth < 0) {
     LOG(INFO) << "Point was behind the camera!";
   }

   LOG(INFO) << "Homogeneous 3D point: " << homogeneous_point3d
             << " reprojected to the pixel value of " << reprojection_pixel;

Point projection can be a tricky function when considering the camera intrinsics
and extrinsics. Theia takes care of this projection (including radial
distortion) in a simple and efficient manner.

In addition to typical getter/setter methods for the camera parameters, the
:class:`Camera` class also defines several helper functions:.

.. function:: bool InitializeFromProjectionMatrix(const int image_width, const int image_height, const Matrix3x4d projection_matrix)

    Initializes the camera intrinsic and extrinsic parameters from the
    projection matrix by decomposing the matrix with a RQ decomposition.

    .. NOTE:: The projection matrix does not contain information about radial
        distortion, so those parameters will need to be set separately.

.. function:: void GetProjectionMatrix(Matrix3x4d* pmatrix) const

    Returns the projection matrix. Does not include radial distortion.

.. function:: void GetCalibrationMatrix(Eigen::Matrix3d* kmatrix) const

    Returns the calibration matrix in the form specified above.

.. function:: Eigen::Vector3d PixelToUnitDepthRay(const Eigen::Vector2d& pixel) const

    Converts the pixel point to a ray in 3D space such that the origin of the
    ray is at the camera center and the direction is the pixel direction rotated
    according to the camera orientation in 3D space. The returned vector is not
    unit length.

Global SfM Pipeline
===================

The global SfM pipelines in Theia follow a general procedure of filtering
outliers and estimating camera poses or structure. Removing outliers can help
increase performance dramatically for global SfM, though robust estimation
methods are still required to obtain good results. The general pipeline is as
follows:

  #. Create the intial view graph from 2-view matches and :class:`TwoViewInfo`
     that describes the relative pose between matched images.
  #. Filter initial view graph and remove outlier 2-view matches.
  #. Calibrate internal parameters of all cameras (either from EXIF or another
     calibration method).
  #. Estimate global orientations of each camera.
  #. Filter the view graph: remove any TwoViewInfos where the relative rotation
     does not agree with the estimated global rotations.
  #. Refine the relative translation estimation to account for the estimated
     global rotations.
  #. Filter any bad :class:`TwoViewInfo` based on the relative translations.
  #. Estimate the global positions of all cameras from the estimated rotations
     and :class:`TwoViewInfo`.
  #. Estimate 3D points.
  #. Bundle adjust the reconstruction.
  #. (Optional) Attempt to estimate any remaining 3D points and bundle adjust again.


Estimating Global Rotations
===========================

Theia estimates the global rotations of cameras robustly using a nonlinear
optimization. Using the relative rotations obtained from all
:class:`TwoViewInfo`, we enforce the constraint that

  :math:`R_{i,j} = R_j * R_i^T`

We use the angle-axis representation of rotations to ensure that proper
rotations are formed. All pairwise constraints are put into a nonlinear
optimization with a robust loss function and the global orienations are
computed. The optimization usually converges within just a few iterations and
provides a very accurate result. The nonlinear optimization is initialized by
forming a random spanning tree of the view graph and walking along the edges.

Estimating Global Positions
===========================

Positions of cameras may be estimated simultaneously after the rotations are
known. We use a nonlinear optimization to estimate camera positions based. Given
pairwise relative translations from :class:`TwoViewInfo` and the estimated
rotation, the constraint

  .. math:: R_i * \dfrac{c_j - c_i}{||c_j - c_i||} = t_{i,j}

Triangulation
=============

  Triangulation in structure from motion calculates the 3D position of an image
  coordinate that has been tracked through several, if not many, images.

  .. cpp:function:: bool Triangulate(const Matrix3x4d& pose1, const Matrix3x4d& pose2, const Eigen::Vector2d& point1, const Eigen::Vector2d& point2, Eigen::Vector4d* triangulated_point)

    2-view triangulation using the method described in [Lindstrom]_. This method
    is optimal in an L2 sense such that the reprojection errors are minimized
    while enforcing the epipolar constraint between the two
    cameras. Additionally, it basically the same speed as the
    :func:`TriangulateDLT` method.

    The poses are the (potentially calibrated) poses of the two cameras, and the
    points are the 2D image points of the matched features that will be used to
    triangulate the 3D point. On successful triangulation, ``true`` is
    returned. The homogeneous 3d point is output so that it may be known if the
    point is at infinity.

  .. function:: bool TriangulateDLT(const Matrix3x4d& pose1, const Matrix3x4d& pose2, const Eigen::Vector2d& point1, const Eigen::Vector2d& point2, Eigen::Vector4d* triangulated_point)

    The DLT triangulation method of [HartleyZisserman]_.

  .. function:: bool TriangulateMidpoint(const Eigen::Vector3d& origin1, const Eigen::Vector3d& ray_direction1, const Eigen::Vector3d& origin2, const Eigen::Vector3d& ray_direction2, Eigen::Vector4d* triangulated_point)

    Perform triangulation by determining the closest point between the two
    rays. In this case, the ray origins are the camera positions and the
    directions are the (unit-norm) ray directions of the features in 3D
    space. This method is known to be suboptimal at minimizing the reprojection
    error, but is approximately 10x faster than the other 2-view triangulation
    methods.

  .. function:: bool TriangulateNViewSVD(const std::vector<Matrix3x4d>& poses, const std::vector<Eigen::Vector2d>& points, Eigen::Vector3d* triangulated_point)
  .. function:: bool TriangulateNView(const std::vector<Matrix3x4d>& poses, const std::vector<Eigen::Vector2d>& points, Eigen::Vector3d* triangulated_point)

    We provide two N-view triangluation methods that minimizes an algebraic
    approximation of the geometric error. The first is the classic SVD method
    presented in [HartleyZisserman]_. The second is a custom algebraic
    minimization. Note that we can derive an algebraic constraint where we note
    that the unit ray of an image observation can be stretched by depth
    :math:`\alpha` to meet the world point :math:`X` for each of the :math:`n`
    observations:

    .. math:: \alpha_i \bar{x_i} = P_i X,

    for images :math:`i=1,\ldots,n`. This equation can be effectively rewritten as:

    .. math:: \alpha_i = \bar{x_i}^\top P_i X,

    which can be substituted into our original constraint such that:

    .. math:: \bar{x_i} \bar{x_i}^\top P_i X = P_i X
    .. math:: 0 = (P_i - \bar{x_i} \bar{x_i}^\top P_i) X

    We can then stack this constraint for each observation, leading to the linear
    least squares problem:

    .. math:: \begin{bmatrix} (P_1 - \bar{x_1} \bar{x_1}^\top P_1) \\ \vdots \\ (P_n - \bar{x_n} \bar{x_n}^\top P_n) \end{bmatrix} X = \textbf{0}

    This system of equations is of the form :math:`AX=0` which can be solved by
    extracting the right nullspace of :math:`A`. The right nullspace of :math:`A`
    can be extracted efficiently by noting that it is equivalent to the nullspace
    of :math:`A^\top A`, which is a 4x4 matrix.

Bundle Adjustment
=================

.. NOTE:: Docmentation coming soon...

Similarity Transformation
=========================

  .. function:: void AlignPointCloudsICP(const int num_points, const double left[], const double right[], double rotation[3 * 3], double translation[3])

    We implement ICP for point clouds. We use Besl-McKay registration to align
    point clouds. We use SVD decomposition to find the rotation, as this is much
    more likely to find the global minimum as compared to traditional ICP, which
    is only guaranteed to find a local minimum. Our goal is to find the
    transformation from the left to the right coordinate system. We assume that
    the left and right reconstructions have the same number of points, and that the
    points are aligned by correspondence (i.e. left[i] corresponds to right[i]).

  .. function:: void AlignPointCloudsUmeyama(const int num_points, const double left[], const double right[], double rotation[3 * 3], double translation[3], double* scale)

    This function estimates the 3D similiarty transformation using the least
    squares method of [Umeyama]_. The returned rotation, translation, and scale
    align the left points to the right such that :math:`Right = s * R * Left +
    t`.

  .. function:: void GdlsSimilarityTransform(const std::vector<Eigen::Vector3d>& ray_origin, const std::vector<Eigen::Vector3d>& ray_direction, const std::vector<Eigen::Vector3d>& world_point, std::vector<Eigen::Quaterniond>* solution_rotation, std::vector<Eigen::Vector3d>* solution_translation, std::vector<double>* solution_scale)

    Computes the solution to the generalized pose and scale problem based on the
    paper "gDLS: A Scalable Solution to the Generalized Pose and Scale Problem"
    by Sweeney et. al. [SweeneyGDLS]_. Given image rays from one coordinate
    system that correspond to 3D points in another coordinate system, this
    function computes the rotation, translation, and scale that will align the
    rays with the 3D points. This is used for applications such as loop closure
    in SLAM and SfM. This method is extremely scalable and highly accurate
    because the cost function that is minimized is independent of the number of
    points. Theoretically, up to 27 solutions may be returned, but in practice
    only 4 real solutions arise and in almost all cases where n >= 6 there is
    only one solution which places the observed points in front of the
    camera. The rotation, translation, and scale are defined such that:
    :math:`sp_i + \alpha_i d_i = RX_i + t` where the observed image ray has an
    origin at :math:`p_i` in the unit direction :math:`d_i` corresponding to 3D
    point :math:`X_i`.

    ``ray_origin``: the origin (i.e., camera center) of the image ray used in
    the 2D-3D correspondence.

    ``ray_direction``: Normalized image rays corresponding to reconstruction
    points. Must contain at least 4 points.

    ``world_point``: 3D location of features. Must correspond to the image_ray
    of the same index. Must contain the same number of points as image_ray, and
    at least 4.

    ``solution_rotation``: the rotation quaternion of the candidate solutions

    ``solution_translation``: the translation of the candidate solutions

    ``solution_scale``: the scale of the candidate solutions
