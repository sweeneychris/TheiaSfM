// Copyright (C) 2014 The Regents of the University of California (Regents).
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

#include <Eigen/Core>
#include <glog/logging.h>
#include <gflags/gflags.h>
#include <theia/theia.h>
#include <string>
#include <vector>

#ifdef __APPLE__
#include <OpenGL/OpenGL.h>
#ifdef FREEGLUT
#include <GL/freeglut.h>
#else
#include <GLUT/glut.h>
#endif
#else
#ifdef _WIN32
#include <windows.h>
#endif
#define GL_GLEXT_PROTOTYPES 1
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#endif

DEFINE_string(reconstruction, "", "Reconstruction file to be viewed.");

// Containers for the data.
std::vector<theia::Camera> cameras;
std::vector<Eigen::Vector3d> world_points;
std::vector<int> num_views_for_track;

// Parameters for OpenGL.
int width = 1200;
int height = 800;
int mouse_down_x[3], mouse_down_y[3];
float rot_x = 0.0f, rot_y = 0.0f;
float prev_x, prev_y;
float distance = 1.0;
Eigen::Vector3d origin = Eigen::Vector3d::Zero();
bool mouse_rotates = false, mouse_moves = false;

bool draw_cameras = true;
bool draw_axes = false;
float point_size = 1.0;
float normalized_focal_length = 1.0;
int min_num_views_for_track = 3;

void GetPerspectiveParams(double* aspect_ratio, double* fovy) {
  double focal_length = 800.0;
  *aspect_ratio = static_cast<double>(width) / static_cast<double>(height);
  *fovy = 2 * atan(height / (2.0 * focal_length)) * 180.0 / M_PI;
}

void ChangeSize(int w, int h) {
  // Prevent a divide by zero, when window is too short
  // (you cant make a window of zero width).
  if (h == 0) h = 1;

  // Use the Projection Matrix
  glMatrixMode(GL_PROJECTION);

  // Reset Matrix
  glLoadIdentity();

  // Set the viewport to be the entire window
  double aspect_ratio, fovy;
  GetPerspectiveParams(&aspect_ratio, &fovy);
  glViewport(0, 0, w, h);

  // Set the correct perspective.
  gluPerspective(fovy, aspect_ratio, 0.001f, 100000.0f);

  // Get Back to the Reconstructionview
  glMatrixMode(GL_MODELVIEW);
}

void DrawAxes(float length) {
  glPushAttrib(GL_POLYGON_BIT | GL_ENABLE_BIT | GL_COLOR_BUFFER_BIT);

  glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
  glDisable(GL_LIGHTING);
  glLineWidth(5.0);
  glBegin(GL_LINES);
  glColor3f(1.0, 0.0, 0.0);
  glVertex3f(0, 0, 0);
  glVertex3f(length, 0, 0);

  glColor3f(0.0, 1.0, 0.0);
  glVertex3f(0, 0, 0);
  glVertex3f(0, length, 0);

  glColor3f(0.0, 0.0, 1.0);
  glVertex3f(0, 0, 0);
  glVertex3f(0, 0, length);
  glEnd();

  glPopAttrib();
  glLineWidth(1.0);
}

void DrawCamera(const theia::Camera& camera) {
  glPushMatrix();
  Eigen::Matrix4d transformation_matrix;
  transformation_matrix.block<3, 3>(0, 0) =
      camera.GetOrientationAsRotationMatrix().transpose();
  transformation_matrix.col(3).head<3>() = camera.GetPosition();
  transformation_matrix(3, 3) = 1.0;

  // Apply world pose transformation.
  glMultMatrixd(reinterpret_cast<GLdouble*>(transformation_matrix.data()));

  // Draw Cameras.
  glColor3f(1.0, 0.0, 0.0);

  // Create the camera wireframe. If intrinsic parameters are not set then use
  // the focal length as a guess.

  double w2f = camera.ImageWidth()/camera.FocalLength()/2.0;
  double h2f = camera.ImageHeight()/camera.FocalLength()/2.0;

  const Eigen::Vector3d top_left =
      normalized_focal_length * Eigen::Vector3d(-w2f, -h2f, 1);
  const Eigen::Vector3d top_right =
      normalized_focal_length * Eigen::Vector3d(w2f, -h2f, 1);
  const Eigen::Vector3d bottom_right =
      normalized_focal_length * Eigen::Vector3d(w2f, h2f, 1);
  const Eigen::Vector3d bottom_left =
      normalized_focal_length * Eigen::Vector3d(-w2f, h2f, 1);

  glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
  glBegin(GL_TRIANGLE_FAN);
  glVertex3f(0.0, 0.0, 0.0);
  glVertex3f(top_right[0], top_right[1], top_right[2]);
  glVertex3f(top_left[0], top_left[1], top_left[2]);
  glVertex3f(bottom_left[0], bottom_left[1], bottom_left[2]);
  glVertex3f(bottom_right[0], bottom_right[1], bottom_right[2]);
  glVertex3f(top_right[0], top_right[1], top_right[2]);
  glEnd();
  glPopMatrix();
}

void RenderScene() {
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glTranslatef(0.0, 0.0, -10);

  glRotatef(180.0f + rot_x, 1.0f, 0.0f, 0.0f);
  glRotatef(-rot_y, 0.0f, 1.0f, 0.0f);
  if (draw_axes) {
    DrawAxes(1000.0);
  }

  glScalef(distance, distance, distance);

  glTranslatef(-origin[0], -origin[1], -origin[2]);
  glClearColor(1.0f, 1.0f, 1.0f, 0.0f);

  // Plot the point cloud.
  glDisable(GL_LIGHTING);
  glEnable(GL_MULTISAMPLE);
  glEnable(GL_BLEND);
  glEnable(GL_POINT_SMOOTH);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  glPointSize(point_size);

  // the coordinates for calculating point attenuation:
  GLfloat point_size_coords[3];
  point_size_coords[0] = 1.0f;
  point_size_coords[1] = 0.055f;
  point_size_coords[2] = 0.0f;
  glPointParameterfv(GL_POINT_DISTANCE_ATTENUATION, point_size_coords);

  // Draw the points.
  glColor3f(0.01, 0.01, 0.01);
  glBegin(GL_POINTS);
  for (int i = 0; i < world_points.size(); i++) {
    if (num_views_for_track[i] < min_num_views_for_track) {
      continue;
    }
    glVertex3d(world_points[i].x(), world_points[i].y(), world_points[i].z());
  }
  glEnd();

  // Draw the cameras.
  if (draw_cameras) {
    for (int i = 0; i < cameras.size(); i++) {
      DrawCamera(cameras[i]);
    }
  }

  glutSwapBuffers();
}

void MouseButton(int button, int state, int x, int y) {
  // button down: save coordinates
  if (state == GLUT_DOWN && button <= 2) {
    mouse_down_x[button] = x;
    mouse_down_y[button] = y;
    prev_x = x;
    prev_y = y;

    if (button == GLUT_RIGHT_BUTTON) {
      mouse_rotates = true;
    } else if (button == GLUT_LEFT_BUTTON) {
      mouse_moves = true;
    }
    return;
  }

  if (state == GLUT_UP && button == GLUT_RIGHT_BUTTON) {
    mouse_rotates = false;
  }

  if (state == GLUT_UP && button == GLUT_LEFT_BUTTON) {
    mouse_moves = false;
  }

  // scroll event - wheel reports as button 3 (scroll up) and button 4 (scroll
  // down)
  if ((button == 3) || (button == 4)) {
    // Each wheel event reports like a button click, GLUT_DOWN then GLUT_UP
    if (state == GLUT_UP) return;  // Disregard redundant GLUT_UP events
    if (button == 3)
      distance /= 1.5f;
    else
      distance *= 1.5f;
  }
}

void MouseMove(int x, int y) {
  if (mouse_rotates) {
    const double rotate_factor = 0.5f;
    // notice x & y difference (i.e., changes in x are to rotate about y-axis)
    rot_x -= rotate_factor * (y - prev_y);
    rot_y -= rotate_factor * (x - prev_x);
    prev_x = x;
    prev_y = y;
  } else if (mouse_moves) {
    const Eigen::Quaterniond inv_rot =
        Eigen::Quaterniond(Eigen::AngleAxisd(theia::DegToRad(180.f + rot_x),
                                             Eigen::Vector3d::UnitX()) *
                           Eigen::AngleAxisd(theia::DegToRad(-rot_y),
                                             Eigen::Vector3d::UnitY()));
    origin += inv_rot * Eigen::Vector3d(prev_x - x, y - prev_y, 0);
    prev_x = x;
    prev_y = y;
  }
}

void Keyboard(unsigned char key, int x, int y) {
  switch (key) {
    case 'r':  // reset viewpoint
      distance = 100.0f;
      rot_x = 0.0f;
      rot_y = 0.0f;
      point_size = 1.0;
      origin = Eigen::Vector3d::Zero();
      break;
    case 'z':
      distance /= 1.2f;
      break;
    case 'Z':
      distance *= 1.2f;
      break;
    case 'p':
      point_size /= 1.2;
      break;
    case 'P':
      point_size *= 1.2;
      break;
    case 'f':
      normalized_focal_length /= 1.2;
      break;
    case 'F':
      normalized_focal_length *= 1.2;
      break;
    case 'c':
      draw_cameras = !draw_cameras;
      break;
    case 'a':
      draw_axes = !draw_axes;
      break;
    case 't':
      ++min_num_views_for_track;
      break;
    case 'T':
      --min_num_views_for_track;
      break;
  }
}

void CenterReconstruction() {
  Eigen::Vector3d mean_camera = Eigen::Vector3d::Zero();
  for (int i = 0; i < cameras.size(); i++) {
    mean_camera += cameras[i].GetPosition();
  }
  mean_camera /= static_cast<double>(cameras.size());

  for (int i = 0; i < cameras.size(); i++) {
    const Eigen::Vector3d old_position = cameras[i].GetPosition();
    cameras[i].SetPosition(old_position - mean_camera);
  }
  for (int i = 0; i < world_points.size(); i++) {
    world_points[i] -= mean_camera;
  }
}

int main(int argc, char* argv[]) {
  THEIA_GFLAGS_NAMESPACE::ParseCommandLineFlags(&argc, &argv, true);
  google::InitGoogleLogging(argv[0]);

  // Output as a binary file.
  std::unique_ptr<theia::Reconstruction> reconstruction(
      new theia::Reconstruction());
  CHECK(ReadReconstruction(FLAGS_reconstruction, reconstruction.get()))
      << "Could not read reconstruction file.";

  // Set up camera drawing.
  cameras.reserve(reconstruction->NumViews());
  for (const theia::ViewId view_id : reconstruction->ViewIds()) {
    const auto* view = reconstruction->View(view_id);
    if (view == nullptr || !view->IsEstimated()) {
      continue;
    }
    cameras.emplace_back(view->Camera());
  }

  // Set up world points.
  world_points.reserve(reconstruction->NumTracks());
  for (const theia::TrackId track_id : reconstruction->TrackIds()) {
    const auto* track = reconstruction->Track(track_id);
    if (track == nullptr || !track->IsEstimated()) {
      continue;
    }
    world_points.emplace_back(track->Point().hnormalized());
    num_views_for_track.emplace_back(track->NumViews());
  }

  CenterReconstruction();

  reconstruction.release();

  // Set up opengl and glut.
  glutInit(&argc, argv);
  glutInitWindowPosition(600, 600);
  glutInitWindowSize(1200, 800);
  glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
  glutCreateWindow("Theia Reconstruction Viewer");

  // Set the camera
  gluLookAt(0.0f, 0.0f, -6.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f);

  // register callbacks
  glutDisplayFunc(RenderScene);
  glutReshapeFunc(ChangeSize);
  glutMouseFunc(MouseButton);
  glutMotionFunc(MouseMove);
  glutKeyboardFunc(Keyboard);
  glutIdleFunc(RenderScene);

  // enter GLUT event processing loop
  glutMainLoop();

  return 0;
}
