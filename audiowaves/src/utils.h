#ifndef _UTILS_H
#define _UTILS_H

#include <cassert>

// ======================================================================

#define square(x) ((x)*(x))
// helper for VBOs
#define BUFFER_OFFSET(i) ((char *)NULL + (i))

#if defined(_WIN32) 
// windows already has them defined...
#define my_max max
#define my_min min
#else
#define my_max std::max
#define my_min std::min
#endif

// ======================================================================

inline glm::vec3 computeNormal(const glm::vec3 &p1, const glm::vec3 &p2, const glm::vec3 &p3) {
  glm::vec3 v12 = p2;
  v12 -= p1;
  glm::vec3 v23 = p3;
  v23 -= p2;
  glm::vec3 normal = glm::cross(v12,v23);
  normal = glm::normalize(normal);
  return normal;
}


inline double triInterpolate(double x_frac, double y_frac, double z_frac,
			     double a, double b, double c, double d, double e, double f, double g, double h) {
  
  assert (x_frac >= 0 && x_frac <= 1);
  assert (y_frac >= 0 && y_frac <= 1);
  assert (z_frac >= 0 && z_frac <= 1);

  // trilinear interpolation
  double ab = (1-z_frac)*a + z_frac*b;
  double cd = (1-z_frac)*c + z_frac*d;
  double ef = (1-z_frac)*e + z_frac*f;
  double gh = (1-z_frac)*g + z_frac*h;
  double abcd = (1-y_frac)*ab + y_frac*cd;
  double efgh = (1-y_frac)*ef + y_frac*gh;
  double abcdefgh = (1-x_frac)*abcd + x_frac*efgh;
  return abcdefgh;
}


inline double AreaOfTriangle(double a, double b, double c) {
  // Area of Triangle =  (using Heron's Formula)
  //  sqrt[s*(s-a)*(s-b)*(s-c)]
  //    where s = (a+b+c)/2
  // also... Area of Triangle = 0.5 * x * c
  double s = (a+b+c) / (double)2;
  double tmp = s*(s-a)*(s-b)*(s-c);
  if (tmp < 0) return 0;
  double answer = sqrt(tmp);
  return answer;
}


inline double AreaOfTriangle(const glm::vec3 &a, const glm::vec3 &b, const glm::vec3 &c) {
  return AreaOfTriangle(glm::length(a-b), glm::length(b-c), glm::length(c-a));
}


void addEdgeGeometry(std::vector<VBOPosNormalColor> &verts,
                     std::vector<VBOIndexedTri> &tri_indices,
                     const glm::vec3 &a, const glm::vec3 &b, 
                     const glm::vec3 &acolor, const glm::vec3 &bcolor, 
                     float a_th,float b_th);



// ======================================================================

#endif

