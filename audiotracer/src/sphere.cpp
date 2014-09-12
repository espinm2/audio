#include "utils.h"
#include "material.h"
#include "argparser.h"
#include "sphere.h"
#include "vertex.h"
#include "mesh.h"
#include "ray.h"
#include "hit.h"

// ====================================================================
// ====================================================================

bool Sphere::intersect(const Ray &r, Hit &h) const {




  // ==========================================
  // ASSIGNMENT:  IMPLEMENT SPHERE INTERSECTION
  // ==========================================

  // plug the explicit ray equation into the implict sphere equation and solve



  // return true if the sphere was intersected, and update the hit
  // data structure to contain the value of t for the ray at the
  // intersection point, the material, and the normal

  // BEGIN SOLUTION
#if 0
  // END SOLUTION
  return false;
  // BEGIN SOLUTION
#endif
  // END SOLUTION


  // BEGIN SOLUTION
  glm::vec3 tmp  = center - r.getOrigin();
  glm::vec3 dir  = r.getDirection();
  float A = glm::dot(dir,dir);
  float B = - 2 * glm::dot(dir,tmp);
  float C = glm::dot(tmp,tmp) - radius*radius;
  float radical = B*B - 4*A*C;
  if (radical < EPSILON) return false;
  
  radical = sqrt(radical);
  float t_m = (-B - radical) / (2*A);
  float t_p = (-B + radical) / (2*A);
  glm::vec3 pt_m = r.pointAtParameter(t_m);
  glm::vec3 pt_p = r.pointAtParameter(t_p);
  
  assert (glm::length(r.getDirection())  > 0.9);

  assert (t_m <= t_p);

  float t = (t_m < EPSILON) ? t_p : t_m;

  if (h.getT() > t  && t > EPSILON) {
    glm::vec3 normal = r.pointAtParameter(t);
    normal -= center;
    normal = glm::normalize(normal);
    h.set(t,this->getMaterial(),normal);
    return true;
  }
  
  return false;
  // END SOLUTION
} 

// ====================================================================
// ====================================================================

// helper function to place a grid of points on the sphere
glm::vec3 ComputeSpherePoint(float s, float t, const glm::vec3 center, float radius) {
  float angle = 2*M_PI*s;
  float y = -cos(M_PI*t);
  float factor = sqrt(1-y*y);
  float x = factor*cos(angle);
  float z = factor*-sin(angle);
  glm::vec3 answer = glm::vec3(x,y,z);
  answer *= radius;
  answer += center;
  return answer;
}

void Sphere::addRasterizedFaces(Mesh *m, ArgParser *args) {
  
  // and convert it into quad patches for radiosity
  int h = args->sphere_horiz;
  int v = args->sphere_vert;
  assert (h % 2 == 0);
  int i,j;
  int va,vb,vc,vd;
  Vertex *a,*b,*c,*d;
  int offset = m->numVertices(); //vertices.size();

  // place vertices
  m->addVertex(center+radius*glm::vec3(0,-1,0));  // bottom
  for (j = 1; j < v; j++) {  // middle
    for (i = 0; i < h; i++) {
      float s = i / float(h);
      float t = j / float(v);
      m->addVertex(ComputeSpherePoint(s,t,center,radius));
    }
  }
  m->addVertex(center+radius*glm::vec3(0,1,0));  // top

  // the middle patches
  for (j = 1; j < v-1; j++) {
    for (i = 0; i < h; i++) {
      va = 1 +  i      + h*(j-1);
      vb = 1 + (i+1)%h + h*(j-1);
      vc = 1 +  i      + h*(j);
      vd = 1 + (i+1)%h + h*(j);
      a = m->getVertex(offset + va);
      b = m->getVertex(offset + vb);
      c = m->getVertex(offset + vc);
      d = m->getVertex(offset + vd);
      m->addRasterizedPrimitiveFace(a,b,d,c,material);
    }
  }

  for (i = 0; i < h; i+=2) {
    // the bottom patches
    va = 0;
    vb = 1 +  i;
    vc = 1 + (i+1)%h;
    vd = 1 + (i+2)%h;
    a = m->getVertex(offset + va);
    b = m->getVertex(offset + vb);
    c = m->getVertex(offset + vc);
    d = m->getVertex(offset + vd);
    m->addRasterizedPrimitiveFace(d,c,b,a,material);
    // the top patches
    va = 1 + h*(v-1);
    vb = 1 +  i      + h*(v-2);
    vc = 1 + (i+1)%h + h*(v-2);
    vd = 1 + (i+2)%h + h*(v-2);
    a = m->getVertex(offset + va);
    b = m->getVertex(offset + vb);
    c = m->getVertex(offset + vc);
    d = m->getVertex(offset + vd);
    m->addRasterizedPrimitiveFace(b,c,d,a,material);
  }
}
