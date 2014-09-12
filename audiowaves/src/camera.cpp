#include "glCanvas.h"
#include "camera.h"

#define _USE_MATH_DEFINES 
#include <cmath>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/matrix_transform.hpp> 
#include <glm/gtx/transform.hpp>
#include <glm/gtx/constants.hpp>

// ====================================================================
// CONSTRUCTORS
// ====================================================================

Camera::Camera(const glm::vec3 &c, const glm::vec3 &poi, const glm::vec3 &u) {
  camera_position = c;
  point_of_interest = poi;
  up = glm::normalize(u);
}

OrthographicCamera::OrthographicCamera
(const glm::vec3 &c, const glm::vec3 &poi, const glm::vec3 &u, double s) 
  : Camera(c,poi,u) {
  size = s;
}

PerspectiveCamera::PerspectiveCamera
(const glm::vec3 &c, const glm::vec3 &poi, const glm::vec3 &u, double a) 
  : Camera(c,poi,u) {
  angle = a;
}

// ====================================================================
// Construct the ViewMatrix & ProjectionMatrix for GL Rendering
// ====================================================================

void OrthographicCamera::glPlaceCamera() {
  glfwGetWindowSize(GLCanvas::window, &width, &height);
  float aspect = width / (float)height;
  float w;
  float h;
  // handle non square windows
  if (aspect < 1.0) {
    w = size / 2.0;
    h = w / aspect;
  } else {
    h = size / 2.0;
    w = h * aspect;
  }
  ProjectionMatrix = glm::ortho<float>(-w,w,-h,h, 0.1f, 100.0f) ;
  ViewMatrix =  glm::lookAt(camera_position,point_of_interest,getScreenUp()) ;
}

void PerspectiveCamera::glPlaceCamera() {
  glfwGetWindowSize(GLCanvas::window, &width, &height);
  float aspect = width / (float)height;
  ProjectionMatrix = glm::perspective<float>(angle, aspect, 0.1f, 100.0f);
  ViewMatrix =  glm::lookAt(camera_position,point_of_interest,getScreenUp()) ;
}

// ====================================================================
// dollyCamera: Move camera along the direction vector
// ====================================================================

void Camera::dollyCamera(double dist) {
  glm::vec3 diff = camera_position - point_of_interest;
  float d = glm::length(diff);
  glm::vec3 translate = float(0.005*d*dist)*getDirection();
  camera_position += translate;
}

// ====================================================================
// zoomCamera: Change the field of view/angle
// ====================================================================

void OrthographicCamera::zoomCamera(double factor) {
  size *= pow(1.005,factor);
}

void PerspectiveCamera::zoomCamera(double dist) {
  angle *= pow(1.003,dist);
}

// ====================================================================
// truckCamera: Translate camera perpendicular to the direction vector
// ====================================================================

void Camera::truckCamera(double dx, double dy) {
  glm::vec3 diff = camera_position - point_of_interest;
  float d = glm::length(diff);
  glm::vec3 translate = (d*0.0007f)*(getHorizontal()*float(dx) + getScreenUp()*float(dy));
  camera_position += translate;
  point_of_interest += translate;
}

// ====================================================================
// rotateCamera: Rotate around the up and horizontal vectors
// ====================================================================

void Camera::rotateCamera(double rx, double ry) {

  // this version of rotate doesn't let the model flip "upside-down"

  // slow the mouse down a little
  rx *= 0.4;
  ry *= 0.4;

  // Note: There is a singularity at the poles (0 & 180 degrees) when
  // 'up' and 'direction' are aligned
  double tiltAngle = acos(glm::dot(up,getDirection())) * 180 / glm::pi<double>();
  if (tiltAngle-ry > 178.0) 
    ry = tiltAngle - 178.0; 
  else if (tiltAngle-ry < 2.0) 
    ry = tiltAngle - 5; 

  glm::vec3 h = getHorizontal();
  glm::mat4 m; 
  m = glm::translate<GLfloat>(m,glm::vec3(point_of_interest));
  m *= glm::rotate<GLfloat>(rx,up);
  m *= glm::rotate<GLfloat>(ry,h);
  m = glm::translate<GLfloat>(m,glm::vec3(-point_of_interest));
  glm::vec4 tmp(camera_position,1);
  tmp = m * tmp;
  camera_position = glm::vec3(tmp.x,tmp.y,tmp.z);
}

// ====================================================================
// ====================================================================

