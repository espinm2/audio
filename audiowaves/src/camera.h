#ifndef _CAMERA_H_
#define _CAMERA_H_

#include <cassert>
#include <iostream>
#include <fstream>

// Include GLM
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

// ====================================================================
// ====================================================================

class Camera {
public:
  // CONSTRUCTOR & DESTRUCTOR
  Camera(const glm::vec3 &c, const glm::vec3 &poi, const glm::vec3 &u);
  virtual ~Camera() {}

  // GL NAVIGATION
  virtual void glPlaceCamera() = 0;
  void dollyCamera(double dist);
  virtual void zoomCamera(double dist) = 0;
  void truckCamera(double dx, double dy);
  void rotateCamera(double rx, double ry);

  const glm::mat4& getViewMatrix() const { return ViewMatrix; }
  const glm::mat4& getProjectionMatrix() const { return ProjectionMatrix; }

public:
  //protected:
  Camera() { assert(0); } // don't use

  // HELPER FUNCTIONS
  glm::vec3 getHorizontal() const {
    return glm::normalize(glm::cross(getDirection(),up));
  }
  glm::vec3 getScreenUp() const {
    return glm::normalize(glm::cross(getHorizontal(),getDirection()));
  }
  glm::vec3 getDirection() const {
    return glm::normalize(point_of_interest - camera_position);
  }

  // REPRESENTATION
  glm::vec3 point_of_interest;
  glm::vec3 camera_position;
  glm::vec3 up;
  int width;
  int height;
  glm::mat4 ViewMatrix;
  glm::mat4 ProjectionMatrix;
};

// ====================================================================

class OrthographicCamera : public Camera {
public:
  // CONSTRUCTOR & DESTRUCTOR
  OrthographicCamera(const glm::vec3 &c = glm::vec3(0,0,1), 
		     const glm::vec3 &poi = glm::vec3(0,0,0), 
		     const glm::vec3 &u = glm::vec3(0,1,0),
		     double s=100);  
  // GL NAVIGATION
  void glPlaceCamera();
  void zoomCamera(double factor);
private:
  // REPRESENTATION
  double size;
};

// ====================================================================

class PerspectiveCamera : public Camera {
public:
  // CONSTRUCTOR & DESTRUCTOR
  PerspectiveCamera(const glm::vec3 &c = glm::vec3(0,0,1), 
		    const glm::vec3 &poi = glm::vec3(0,0,0), 
		    const glm::vec3 &u = glm::vec3(0,1,0),
		    double a = 45);
  // GL NAVIGATION
  void glPlaceCamera();
  void zoomCamera(double dist);
private:
  // REPRESENTATION
  double angle;
};

// ====================================================================

#endif

