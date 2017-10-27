// Copyright (C) 2016 The Regents of the University of California (Regents).
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//
//     * Neither the name of The Regents or University of California nor the
//       names of its contributors may be used to endorse or promote products
//       derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Please contact the author of this library if you have any questions.
// Author: Chris Sweeney (cmsweeney@cs.ucsb.edu)

#include "theia/sfm/undistort_image.h"

#include <Eigen/Core>
#include <algorithm>
#include <limits>

#include "theia/image/image.h"
#include "theia/sfm/camera/camera.h"
#include "theia/sfm/camera/camera_intrinsics_model.h"
#include "theia/sfm/camera/division_undistortion_camera_model.h"
#include "theia/sfm/camera/fisheye_camera_model.h"
#include "theia/sfm/camera/fov_camera_model.h"
#include "theia/sfm/camera/pinhole_camera_model.h"
#include "theia/sfm/camera/pinhole_radial_tangential_camera_model.h"
#include "theia/sfm/reconstruction.h"
#include "theia/sfm/track.h"
#include "theia/sfm/types.h"
#include "theia/sfm/view.h"

namespace theia {
namespace {

// Set all radial distortion parameters to zero. This requires adjusting
// different parameters for each camera intrinsics model type.
void SetLensDistortionToZero(Camera* camera) {
  double* intrinsics = camera->mutable_intrinsics();
  switch (camera->GetCameraIntrinsicsModelType()) {
    case CameraIntrinsicsModelType::PINHOLE:
      intrinsics[PinholeCameraModel::RADIAL_DISTORTION_1] = 0.0;
      intrinsics[PinholeCameraModel::RADIAL_DISTORTION_2] = 0.0;
      break;
    case CameraIntrinsicsModelType::PINHOLE_RADIAL_TANGENTIAL:
      intrinsics[PinholeRadialTangentialCameraModel::RADIAL_DISTORTION_1] = 0.0;
      intrinsics[PinholeRadialTangentialCameraModel::RADIAL_DISTORTION_2] = 0.0;
      intrinsics[PinholeRadialTangentialCameraModel::RADIAL_DISTORTION_3] = 0.0;
      intrinsics[PinholeRadialTangentialCameraModel::TANGENTIAL_DISTORTION_1] =
          0.0;
      intrinsics[PinholeRadialTangentialCameraModel::TANGENTIAL_DISTORTION_2] =
          0.0;
      break;
    case CameraIntrinsicsModelType::FISHEYE:
      intrinsics[FisheyeCameraModel::RADIAL_DISTORTION_1] = 0.0;
      intrinsics[FisheyeCameraModel::RADIAL_DISTORTION_2] = 0.0;
      intrinsics[FisheyeCameraModel::RADIAL_DISTORTION_3] = 0.0;
      intrinsics[FisheyeCameraModel::RADIAL_DISTORTION_4] = 0.0;
      break;
    case CameraIntrinsicsModelType::FOV:
      intrinsics[FOVCameraModel::RADIAL_DISTORTION_1] = 0.0;
      break;
    case CameraIntrinsicsModelType::DIVISION_UNDISTORTION:
      intrinsics[DivisionUndistortionCameraModel::RADIAL_DISTORTION_1] = 0.0;
      break;
    default:
      LOG(FATAL) << "Invalid camera intrinsics type.";
      break;
  }
}

// We seek a mapping from undistorted pixels to distorted pixels such that a
// valid mapping exists for all undistorted pixels. This is effectively
// undistorting the image then cropping it so that no unmapped (black) pixels
// are in the undistorted image. To do this, we need to find the min and max x/y
// values of unmapped pixels so that we can determine the size and location of
// the crop.
void FindUndistortedImageBoundary(const Camera& distorted_camera,
                                  const Camera& undistorted_camera,
                                  Eigen::Vector4d* bounds) {
  const CameraIntrinsicsModel& distorted_intrinsics =
      *distorted_camera.CameraIntrinsics();
  const CameraIntrinsicsModel& undistorted_intrinsics =
      *undistorted_camera.CameraIntrinsics();

  // Find the max and min locations of the undistorted pixels.
  double left_max_x = std::numeric_limits<double>::lowest();
  double right_min_x = std::numeric_limits<double>::max();
  for (size_t y = 0; y < distorted_camera.ImageHeight(); ++y) {
    // Left border.
    const Eigen::Vector3d distorted_point1 =
        distorted_intrinsics.ImageToCameraCoordinates(
            Eigen::Vector2d(0.5, y + 0.0));
    const Eigen::Vector2d undistorted_point1 =
        undistorted_intrinsics.CameraToImageCoordinates(distorted_point1);
    left_max_x = std::max(left_max_x, undistorted_point1(0));

    // Right border.
    const Eigen::Vector3d distorted_point2 =
        distorted_intrinsics.ImageToCameraCoordinates(
            Eigen::Vector2d(distorted_camera.ImageWidth() - 0.5, y + 0.5));
    const Eigen::Vector2d undistorted_point2 =
        undistorted_intrinsics.CameraToImageCoordinates(distorted_point2);
    right_min_x = std::min(right_min_x, undistorted_point2(0));
  }

  // Determine min, max coordinates along left / right image border.
  double top_max_y = std::numeric_limits<double>::lowest();
  double bottom_min_y = std::numeric_limits<double>::max();
  for (size_t x = 0; x < distorted_camera.ImageWidth(); ++x) {
    // Top border.
    const Eigen::Vector3d distorted_point1 =
        distorted_intrinsics.ImageToCameraCoordinates(
            Eigen::Vector2d(x + 0.5, 0.5));
    const Eigen::Vector2d undistorted_point1 =
        undistorted_intrinsics.CameraToImageCoordinates(distorted_point1);
    top_max_y = std::max(top_max_y, undistorted_point1(1));

    // Bottom border.
    const Eigen::Vector3d distorted_point2 =
        distorted_intrinsics.ImageToCameraCoordinates(
            Eigen::Vector2d(x + 0.5, distorted_camera.ImageHeight() - 0.5));
    const Eigen::Vector2d undistorted_point2 =
        undistorted_intrinsics.CameraToImageCoordinates(distorted_point2);
    bottom_min_y = std::min(bottom_min_y, undistorted_point2(1));
  }

  *bounds = Eigen::Vector4d(left_max_x, right_min_x, top_max_y, bottom_min_y);
}

// Create an undistorted image from the distorted image given the distorted and
// undistorted camera parameters. This function only maps the pixels to create
// the undistorted image and assumes the camera parameters and image sizes have
// already been solved for.
void RemoveImageLensDistortion(const Camera& distorted_camera,
                               const FloatImage& distorted_image,
                               const Camera& undistorted_camera,
                               FloatImage* undistorted_image) {
  const CameraIntrinsicsModel& distorted_intrinsics =
      *distorted_camera.CameraIntrinsics();
  const CameraIntrinsicsModel& undistorted_intrinsics =
      *undistorted_camera.CameraIntrinsics();

  // For each pixel in the undistorted image, find the coordinate in the
  // distorted image and set the pixel color accordingly.
  const int num_channels = distorted_image.Channels();
  OpenImageIO::ImageBuf& undistorted_img =
      undistorted_image->GetOpenImageIOImageBuf();
  OpenImageIO::ImageBuf::Iterator<float> undistorted_it(undistorted_img);
  for (; !undistorted_it.done(); ++undistorted_it) {
    Eigen::Vector2d image_point(undistorted_it.x() + 0.5,
                                undistorted_it.y() + 0.5);

    // Camera models assume that the upper left pixel center is (0.5, 0.5).
    const Eigen::Vector3d distorted_point =
        undistorted_intrinsics.ImageToCameraCoordinates(image_point);
    const Eigen::Vector2d distorted_pixel =
        distorted_intrinsics.CameraToImageCoordinates(distorted_point);
    const Eigen::Vector2d pixel(std::round<int>(distorted_pixel.x() - 0.5),
                                std::round<int>(distorted_pixel.y() - 0.5));
    // Set all color channels appropriately. Note that we do not need to check
    // if the distorted pixel is within the image bounds since the
    // FindUndistortedImageBoundary function has guaranteed that all pixels will
    // be within the image borders.
    for (int c = 0; c < num_channels; c++) {
      undistorted_it[c] = distorted_image.BilinearInterpolate(
          distorted_pixel.x(), distorted_pixel.y(), c);
    }
  }
}

}  // namespace

bool UndistortImage(const Camera& distorted_camera,
                    const FloatImage& distorted_image,
                    const Camera& undistorted_camera,
                    FloatImage* undistorted_image) {
  // Undistort the image.
  if (distorted_image.Channels() == 1) {
    undistorted_image->ConvertToGrayscaleImage();
  } else {
    undistorted_image->ConvertToRGBImage();
  }
  undistorted_image->Resize(undistorted_camera.ImageWidth(),
                            undistorted_camera.ImageHeight());

  // Remap the distorted pixels into the undistorted image.
  RemoveImageLensDistortion(
      distorted_camera, distorted_image, undistorted_camera, undistorted_image);

  return true;
}

// Create the undistorted camera by removing radial distortion parameters.
bool UndistortCamera(const Camera& distorted_camera,
                     Camera* undistorted_camera) {
  *undistorted_camera = distorted_camera;
  SetLensDistortionToZero(undistorted_camera);

  Eigen::Vector4d undistorted_image_boundaries;
  FindUndistortedImageBoundary(
      distorted_camera, *undistorted_camera, &undistorted_image_boundaries);

  // Given the locations of the min/max undistorted pixels, compute the scale
  // factor to resize the undistorted image.
  const double cx = undistorted_camera->PrincipalPointX();
  const double cy = undistorted_camera->PrincipalPointY();
  const double left_max_x = undistorted_image_boundaries(0);
  const double right_min_x = undistorted_image_boundaries(1);
  const double top_max_y = undistorted_image_boundaries(2);
  const double bottom_min_y = undistorted_image_boundaries(3);

  // Scale undistorted camera dimensions.
  const double scale_x =
      1.0 /
      std::max(cx / (cx - left_max_x),
               (distorted_camera.ImageWidth() - 0.5 - cx) / (right_min_x - cx));
  const double scale_y =
      1.0 / std::max(cy / (cy - top_max_y),
                     (distorted_camera.ImageHeight() - 0.5 - cy) /
                         (bottom_min_y - cy));

  undistorted_camera->SetImageSize(
      static_cast<int>(
          std::max(1.0, scale_x * undistorted_camera->ImageWidth())),
      static_cast<int>(
          std::max(1.0, scale_y * undistorted_camera->ImageHeight())));

  // Scale the principal point according to the new dimensions of the image.
  undistorted_camera->SetPrincipalPoint(
      undistorted_camera->PrincipalPointX() *
          static_cast<double>(undistorted_camera->ImageWidth()) /
          distorted_camera.ImageWidth(),
      undistorted_camera->PrincipalPointY() *
          static_cast<double>(undistorted_camera->ImageHeight()) /
          distorted_camera.ImageHeight());

  return true;
}

bool UndistortReconstruction(Reconstruction* reconstruction) {
  const auto view_ids = reconstruction->ViewIds();
  for (const ViewId view_id : view_ids) {
    View* view = reconstruction->MutableView(view_id);
    if (view == nullptr || !view->IsEstimated()) {
      continue;
    }

    // Undistort the image to obtain the undistorted camera.
    const Camera distorted_camera = view->Camera();
    Camera* undistorted_camera = view->MutableCamera();
    if (!UndistortCamera(distorted_camera, undistorted_camera)) {
      return false;
    }

    // The camera intrinsics models describe how to distort and undistort the
    // pixels.
    const CameraIntrinsicsModel& distorted_intrinsics =
        *distorted_camera.CameraIntrinsics();
    const CameraIntrinsicsModel& undistorted_intrinsics =
        *undistorted_camera->CameraIntrinsics();

    // Undisort all features seen by the view.
    const auto track_ids = view->TrackIds();
    for (const TrackId track_id : track_ids) {
      Track* track = reconstruction->MutableTrack(track_id);
      if (track == nullptr || !track->IsEstimated()) {
        continue;
      }

      // Get the undistorted feature.
      const Feature* distorted_feature = view->GetFeature(track_id);
      const Eigen::Vector3d distorted_point =
          distorted_intrinsics.ImageToCameraCoordinates(*distorted_feature);
      const Eigen::Vector2d undistorted_feature =
          undistorted_intrinsics.CameraToImageCoordinates(distorted_point);

      // Add the new undistorted feature to the undistorted view. First, remove
      // the feature since it is not mutable.
      CHECK(view->RemoveFeature(track_id));

      // Only add the feature back if it is within the new image
      // boundaries. Otherwise, remove the feature from the view and the track.
      if (undistorted_feature.x() < 0 ||
          undistorted_feature.x() > undistorted_camera->ImageWidth() - 1 ||
          undistorted_feature.y() < 0 ||
          undistorted_feature.y() > undistorted_camera->ImageHeight() - 1) {
        track->RemoveView(view_id);
        // Remove the track if there are no more observations.
        if (track->NumViews() == 0) {
          reconstruction->RemoveTrack(track_id);
        }
      } else {
        view->AddFeature(track_id, undistorted_feature);
      }
    }
  }

  return true;
}

}  // namespace theia
