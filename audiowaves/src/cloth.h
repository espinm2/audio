#ifndef _CLOTH_H_
#define _CLOTH_H_

#include "argparser.h"
#include "boundingbox.h"
#include "vbo_structs.h"
#include <vector>
#include <map>

using std::vector;
using std::pair;

// =====================================================================================
// Cloth Particles
// =====================================================================================

class ClothParticle {
public:

  // ACCESSORS
  const Vec3f& getOriginalPosition() const{ return original_position; }
  const Vec3f& getPosition() const{ return position; }
  const Vec3f& getVelocity() const{ return velocity; }
  const Vec3f& getAcceleration() const { return acceleration; }
  Vec3f getForce() const { return mass*acceleration; }
  double getMass() const { return mass; }
  bool isFixed() const { return fixed; }
  bool isLoose() const {return !fixed; }

  // MODIFIERS
  void setOriginalPosition(const Vec3f &p) { original_position = p; }
  void setPosition(const Vec3f &p) { position = p; }
  void setVelocity(const Vec3f &v) { velocity = v; }
  void setAcceleration(const Vec3f &a) { acceleration = a; }
  void setMass(double m) { mass = m; }
  void setFixed(bool b) { fixed = b; }
private:
  // REPRESENTATION
  Vec3f original_position;
  Vec3f position;
  Vec3f velocity;
  Vec3f acceleration;
  double mass;
  bool fixed;

};

// =====================================================================================
// Cloth System
// =====================================================================================


class Cloth {

public:
  Cloth(ArgParser *args);
  ~Cloth() { delete [] particles; cleanupVBOs(); }

  // ACCESSORS
  const BoundingBox& getBoundingBox() const { return box; }
  const vector<ClothParticle*> getAdjParticles (int i, int j) ;
  const vector<ClothParticle*> getShearParticles (int i, int j) ;
  const vector<ClothParticle*> getFlexParticles (int i, int j) ;
  const Vec3f getSpringForce(ClothParticle* a, ClothParticle* b) ;

  // PAINTING & ANIMATING
  void Paint() const; // <------What is this Paint?
  void Animate();

  void initializeVBOs();
  void setupVBOs();
  void drawVBOs();
  void cleanupVBOs();

private:

  // PRIVATE ACCESSORS
  const ClothParticle& getParticle(int i, int j) const {
    assert (i >= 0 && i < nx && j >= 0 && j < ny);
    return particles[i + j*nx]; }

  ClothParticle& getParticle(int i, int j) {
    assert (i >= 0 && i < nx && j >= 0 && j < ny);
    return particles[i + j*nx]; }

  Vec3f computeGouraudNormal(int i, int j) const;

  // HELPER FUNCTION
  void computeBoundingBox();
  void AddVBOEdge(int i1, int j1, int i2, int j2, double correction);

  // REPRESENTATION
  ArgParser *args;
  // grid data structure
  int nx, ny;
  ClothParticle *particles;
  BoundingBox box;
  // simulation parameters
  double damping;
  // spring constants
  double k_structural;
  double k_shear;
  double k_bend;
  // correction thresholds
  double provot_structural_correction;
  double provot_shear_correction;


  // VBOs
  GLuint cloth_verts_VBO;
  GLuint cloth_quad_indices_VBO;
  GLuint cloth_happy_edge_indices_VBO;
  GLuint cloth_unhappy_edge_indices_VBO;
  GLuint cloth_velocity_visualization_VBO;
  GLuint cloth_force_visualization_VBO;
  std::vector<VBOPosNormalColor> cloth_verts; 
  std::vector<VBOIndexedQuad> cloth_quad_indices;
  std::vector<VBOIndexedEdge> cloth_happy_edge_indices;
  std::vector<VBOIndexedEdge> cloth_unhappy_edge_indices;
  std::vector<VBOPosColor> cloth_velocity_visualization;
  std::vector<VBOPosColor> cloth_force_visualization;

};

// ========================================================================

#endif
