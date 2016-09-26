.. highlight:: c++

.. default-domain:: cpp

.. _`chapter-cameras`:

=============
Camera Models
=============

At the base of Structure-from-Motion is the imaging device: a camera. Cameras
capture the pixels that are the input to SfM. Most consumer point-and-shoot
cameras are standard pinhole perspective cameras that may have slight radial
lens distortion; however, other cameras such as fisheye and omnidirectional
cameras exist and may be useful for SfM due to their wide fields of view. In
Theia, we allow for any camera model so long as the projection from 3D space on
to the imaging plane (and the reverse transformation) is well-defined. Theia
uses polymorphism so that new cameras that are added to the library are
seamlessly integrated with the SfM pipelines. This greatly simplifies the SfM
code, and make it simple for users to add new camera models.

First, it is useful to review the coordinate system conventions in Theia. We
utilize 3 coordinate systems: world, camera, and image coordinate systems. The
world coordinate system is the global coordinate system of 3D space that defines
3D point locations and camera positions. The camera coordinate system is
cenetered at the particular camera of interest and is oriented so that the
positive z-axis is aligned with the camera's optical axis (i.e. the optical axis
is the ray :math:`\left[0 & 0 & 1]` in the camera coordinates. The image
coordinate system is the 2D coordinate system that describes image pixel
coordinates. The origin is at the top-left of the image with positive x going
towards the right and positive y going down. All coordinate systems are
right-handed.

At the top level, Theia contains a :class:`Camera` class. This class contains
the camera's extrinsics pose (i.e., the orientation and position in 3D space) as
well as a projection model that defines how the camera projects 3D points onto
the image plane (via the :class:`CameraIntrinsicsModel` class). The type of
projection model is defined at runtime, and the :class:`Camera` API is agnostic
to the projection model. Currently there are 4 types of projection models:
:class:`PinholeCameraModel`, :class:`PinholeRadialTangentialCameraModel`,
:class:`FisheyeCameraModel`, and :class:`FOVCameraModel`.

Camera
------

.. class:: Camera

The :class:`Camera` class contains intrinsic and extrinsic information about the
camera that observed the scene. Theia has an efficient, compact :class:`Camera`
class that abstracts away common image operations. This greatly relieves the
pain of manually dealing with calibration and geometric transformations of
images. The projection model (i.e. perspective, fisheye, etc.) is defined by the
user at runtime.

We store the camera pose information as the transformation which maps world
coordinates into camera coordinates. Our rotation is stored internally as an
angle-axis rotation, which makes optimization with :class:`BundleAdjustment`
more effective. However, for convenience we provide an interface to retrieve the
rotation as a rotation matrix as well. Further, we store the camera position as
opposed to the translation.

The convenience of this camera class is clear with the common example of 3D
point reprojection.

.. code:: c++

   // Open an image and obtain camera parameters.
   FloatImage image("my_image.jpg");
   const Eigen::Matrix3d rotation = value obtained elsewhere...
   const Eigen::Vector3d position = value obtained elsewhere...

   // Set up the camera.
   Camera camera;
   camera.SetOrientationFromRotationMatrix(rotation);
   camera.SetPosition(position);

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
and extrinsics. Theia provides the convenient API for these sorts of functions
that affords users a clean interface and the ability to mix and match various
camera models.

In addition to typical getter/setter methods for the camera parameters, the
:class:`Camera` class also defines several helper functions:.

.. function:: void SetFromCameraIntrinsicsPriors(const CameraIntrinsicsPrior& prior)

    Sets the camera intrinsics parameters from the priors, including the camera model.

.. function:: bool Camera::InitializeFromProjectionMatrix(const int image_width, const int image_height, const Matrix3x4d projection_matrix)

    Initializes the camera intrinsic and extrinsic parameters from the
    projection matrix by decomposing the matrix with a RQ decomposition.

    .. NOTE:: The projection matrix does not contain information about radial
        distortion, so those parameters will need to be set separately.

.. function:: void Camera::GetProjectionMatrix(Matrix3x4d* pmatrix) const

    Returns the projection matrix. Does not include radial distortion.

.. function:: void Camera::GetCalibrationMatrix(Eigen::Matrix3d* kmatrix) const

    Returns the calibration matrix in the form specified above.

.. function:: Eigen::Vector3d Camera::PixelToUnitDepthRay(const Eigen::Vector2d& pixel) const

    Converts the pixel point to a ray in 3D space such that the origin of the
    ray is at the camera center and the direction is the pixel direction rotated
    according to the camera orientation in 3D space. The returned vector is not
    unit length.


CameraIntrinsicsModel
---------------------

.. class:: CameraIntrinsicsModel

The projection of 3D points into image pixels is defined by the camera
model. This model depends on the type of lens being used, the field of view, and
more. Different camera models have different benefits: most consumer cameras may
be modelled with perspective projection, but wide field of view cameras such as
GoPros are modelled more appropriately with a fishey camera model. To allow for
any type of camera projection and distortion, Theia utilizes an abstract
interface :class:`CameraIntrinsicsModel` class. This class defines the interface
for projection and un-projection, as well as several methods other that subclasses are
required to implement.

.. function:: CameraIntrinsicsModelType CameraIntrinsicsModel::Type()

    Each camera intrinsics model that is implemented will have a type (found in the enum :class:`CameraIntrinsicsModelType` in camera_intrinsics_model_type.h. This type is unique to each implemented camera model

.. function:: int CameraIntrinsicsModel::NumParameters()

    Returns the number of camera intrinsics parameters that are used for the
    particular camera model. This is the number of "free" parameters (i.e., ones
    that may be optimized) for the camera model.

.. function:: void CameraIntrinsicsModel::SetFromCameraIntrinsicsPrior()

    The :class:`CameraIntrinsicsPrior` class specifies metadata and prior
    information that may be used to initialize camera parameters. For example,
    this class may contain a focal length extracted from EXIF metadata.

.. function:: CameraIntrinsicsPrior CameraIntrinsicsModel::CameraIntrinsicsPriorFromIntrinsics()

    Returns a CameraIntrinsicsPrior object populated with the appropriate fields
    related to the camera intrinsic parameters.

.. function:: void CameraIntrinsicsModel::GetSubsetFromOptimizeIntrinsicsType(const OptimizeIntrinsicsType& intrinsics_to_optimize)

    :class:`BundleAdjustment` allows for individual camera parameters to be optimized or set constant. Since each derived :class:`CameraIntrinsicsModel` class may contain different intrinsics, this helper method returns the appropriate indices of parameters that should be kept constant during optimization based on the intrinsics_to_optimize input.

.. function:: Eigen::Vector2d CameraIntrinsicsModel::CameraToImageCoordinates(const Eigen::Vector3d& point)

    Projects the 3D point in the camera coordinate system (NOTE: this is
    different from the "world coordinate system") into the image
    coordinates. This includes apply lens/radial distortion.

.. function:: Eigen::Vector3d CameraIntrinsicsModel::ImageToCameraCoordinates(const Eigen::Vector2d& pixel)

    Given the pixel coordinate, this method returns the ray corresponding to the
    pixel. This involves removing the effects of camera intrinsics and
    lens distortion.

.. function:: Eigen::Vector2d CameraIntrinsicsModel::DistortPoint(const Eigen::Vector2d& point)

    Given the point in camera coordinates, apply lens distortion.

.. function:: Eigen::Vector2d CameraIntrinsicsModel::UndistortPoint(const Eigen::Vector2d& point)

    Given the distorted point in camera coordinates, remove the effects of lens distortion.


PinholeCameraModel
---------------------

.. class:: PinholeCameraModel

The Pinhole camera model is the most common camera model for consumer
cameras. In this model, the image is mapped onto a plane through perspective
projection. The projection is defined by the camera intrinsic parameters such as
focal length, principal point, aspect ratio, and skew. These parameters define
an intrinsics matrix:

.. math::
  K = \left[\begin{matrix}f & s & p_x \\ 0 & f * a & p_y \\ 0 & 0 & 1 \end{matrix} \right]

where :math:`f` is the focal length (in pixels), :math:`s` is the skew,
:math:`a` is the aspect ratio and :math:`p` is the principle point of the
camera. All of these intrinsics may be accessed with getter and setter methods,
e.g., :code:`double GetFocalLength()` or :code:`void SetFocalLength(const double
focal_length)`. Note that we do additionally allow for up to two radial
distortion parameters that model lens distortion.

.. class:: PinholeRadialTangentialCameraModel

This class is the same as the :class:`PinholeCameraModel` but includes 3 radial
distortion and 2 tangential distortion parameters.


FisheyeCameraModel
------------------

.. class:: FisheyeCameraModel

The Fisheye camera model is a camera model utilized for wide field of view
cameras. This camera model is neccessary because the pinhole perspective camera
model is not capable of modeling image projections as the field of view
approaches 180 degrees. The camera model is based on the `OpenCV fisheye camera model <http://docs.opencv.org/2.4/modules/calib3d/doc/camera_calibration_and_3d_reconstruction.html#fisheye>`_

Given a point :math:`X=\left[\begin{matrix}x & y & z\end{matrix} \right]` in camera coordinates, the fisheye projection is:

.. math::

    r = \sqrt{x^2 + y^2}
    \theta = atan2(r, |z|) \\
    \theta_d = \theta (1 + k1 * \theta^2 + k2 * \theta^4 + k3 * \theta^6 + k4 * \theta^8) \\
    x' = \theta_d * x / r \\
    y' = \theta_d * y / r \\

Where :math:`\left[x'  y' \right]` is the projected (and distorted) image
point. This projection model uses the angle between the observed point and the
camera's optical axis to determine the projection and the distortion. This allows for observations near or above the 180 degree field of view.

FovCameraModel
--------------

.. class:: FOVCameraModel

This class contains the camera intrinsic information for fov cameras. This is
an alternative representation for camera models with large radial distortion
(such as fisheye cameras) where the distance between an image point and
principal point is roughly proportional to the angle between the 3D point and
the optical axis. This camera model is first proposed in [Devernay]_.


Adding a New Camera Model
-------------------------

The CameraIntrinsicsModel describes the abstract interface for mapping between camera and image
coordinate systems. To implement a new camera model, you will have to take the following steps.

1) Create a derived class from this :class:`CameraIntrinsicsModel`, and
   implement all of the pure virtual methods and the static methods that are
   used for camera projection.

2) Add an enum to :class:`CameraIntrinsicsModelType` and add an "else if" to the
   :func:`Create` method in this class to allow your camera model to be created.

3) Add the new class and its :class:`CameraIntrinsicsType` to the
   CAMERA_MODEL_SWITCH_STATEMENT macro in camera_intrinsics_model.cc

4) Add a switch/case in create_reprojection_error_cost_function.h to handle
   the new camera model.

5) Create unit tests to ensure that your new camera model is functioning
   properly!
