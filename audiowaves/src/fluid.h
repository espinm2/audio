#ifndef _FLUID_H_
#define _FLUID_H_

#include <cassert>
#include <vector>
#include "argparser.h"
#include "boundingbox.h"
#include "cell.h"
#include "vbo_structs.h"

class ArgParser;
class MarchingCubes;

// ========================================================================
// ========================================================================

class Fluid {

public:

  // ========================
  // CONSTRUCTOR & DESTRUCTOR
  Fluid(ArgParser *_args);
  ~Fluid();

  // load data
  void Load();

  // ===============================
  // ANIMATION & RENDERING FUNCTIONS
  void Animate();
  void initializeVBOs(); 
  void setupVBOs(); 
  void drawVBOs();
  void cleanupVBOs();
  BoundingBox getBoundingBox() const {
    return BoundingBox(glm::vec3(0,0,0),glm::vec3(nx*dx,ny*dy,nz*dz)); }

private:

  // ==============
  // CELL ACCESSORS
  int Index(int i, int j, int k) const {
    assert (i >= -1 && i <= nx);
    assert (j >= -1 && j <= ny);
    assert (k >= -1 && k <= nz);
    return (i+1)*(ny+2)*(nz+2) + (j+1)*(nz+2) + (k+1);
  }
  Cell* getCell(int i, int j, int k) const { return &cells[Index(i,j,k)]; }

  // =================
  // ANIMATION HELPERS
  void ComputeNewVelocities();
  void SetBoundaryVelocities();
  void EmptyVelocities(int i, int j, int k);
  void CopyVelocities();
  double AdjustForIncompressibility();
  void UpdatePressures();
  void MoveParticles();
  void ReassignParticles();
  void SetEmptySurfaceFull();

  // =====================
  // NAVIER-STOKES HELPERS
  glm::vec3 getInterpolatedVelocity(const glm::vec3 &pos) const;
  double getPressure(int i, int j, int k) const { return getCell(i,j,k)->getPressure(); }
  // velocity accessors
  double get_u_plus(int i, int j, int k) const { return getCell(i,j,k)->get_u_plus(); }
  double get_v_plus(int i, int j, int k) const { return getCell(i,j,k)->get_v_plus(); }  
  double get_w_plus(int i, int j, int k) const { return getCell(i,j,k)->get_w_plus(); }  
  double get_new_u_plus(int i, int j, int k) const { return getCell(i,j,k)->get_new_u_plus(); }  
  double get_new_v_plus(int i, int j, int k) const { return getCell(i,j,k)->get_new_v_plus(); }  
  double get_new_w_plus(int i, int j, int k) const { return getCell(i,j,k)->get_new_w_plus(); }  
  double get_u_avg(int i, int j, int k) const { return 0.5*(get_u_plus(i-1,j,k)+get_u_plus(i,j,k)); }
  double get_v_avg(int i, int j, int k) const { return 0.5*(get_v_plus(i,j-1,k)+get_v_plus(i,j,k)); }
  double get_w_avg(int i, int j, int k) const { return 0.5*(get_w_plus(i,j,k-1)+get_w_plus(i,j,k)); }
  double get_uv_plus(int i, int j, int k) const { 
    return 0.5*(get_u_plus(i,j,k) + get_u_plus(i,j+1,k)) * 0.5*(get_v_plus(i,j,k) + get_v_plus(i+1,j,k)); }
  double get_uw_plus(int i, int j, int k) const { 
    return 0.5*(get_u_plus(i,j,k) + get_u_plus(i,j,k+1)) * 0.5*(get_w_plus(i,j,k) + get_w_plus(i+1,j,k)); }
  double get_vw_plus(int i, int j, int k) const { 
    return 0.5*(get_v_plus(i,j,k) + get_v_plus(i,j,k+1)) * 0.5*(get_w_plus(i,j,k) + get_w_plus(i,j+1,k)); }
  // velocity modifiers
  void set_new_u_plus(int i, int j, int k, double f) { getCell(i,j,k)->set_new_u_plus(f); }
  void set_new_v_plus(int i, int j, int k, double f) { getCell(i,j,k)->set_new_v_plus(f); }
  void set_new_w_plus(int i, int j, int k, double f) { getCell(i,j,k)->set_new_w_plus(f); }
  void adjust_new_u_plus(int i, int j, int k, double f) { getCell(i,j,k)->adjust_new_u_plus(f); }
  void adjust_new_v_plus(int i, int j, int k, double f) { getCell(i,j,k)->adjust_new_v_plus(f); }
  void adjust_new_w_plus(int i, int j, int k, double f) { getCell(i,j,k)->adjust_new_w_plus(f); }

  // ========================================
  // RENDERING SURFACE (using Marching Cubes)
  double interpolateIsovalue(const glm::vec3 &c) const;
  double getIsovalue(int i, int j, int k) const;

  // ============
  // LOAD HELPERS
  bool inShape(glm::vec3 &pos, const std::string &shape);
  void GenerateParticles(const std::string &shape, const std::string &placement);

  // don't use this constructor
  Fluid() { assert(0); }
  
  // ==============
  // REPRESENTATION
  ArgParser *args;

  // fluid parameters
  int nx,ny,nz;     // number of grid cells in each dimension
  double dx,dy,dz;  // dimensions of each grid cell
  Cell *cells;      // NOTE: padded with extra cells on each side

  // simulation parameters
  bool xy_free_slip;
  bool yz_free_slip;
  bool zx_free_slip;
  bool compressible;
  double viscosity;
  double density; // average # of particles initialized in each "Full" cell

  // VBOs
  GLuint fluid_particles_VBO;
  GLuint fluid_velocity_verts_VBO;
  GLuint fluid_velocity_tri_indices_VBO;
  GLuint fluid_facevelocity_verts_VBO;
  GLuint fluid_facevelocity_tri_indices_VBO;
  GLuint fluid_pressure_vis_VBO;
  GLuint fluid_cell_type_vis_VBO;

  std::vector<VBOPosNormalColor> fluid_particles;
  std::vector<VBOPosNormalColor> fluid_velocity_verts; 
  std::vector<VBOIndexedTri> fluid_velocity_tri_indices;
  std::vector<VBOPosNormalColor> fluid_facevelocity_verts; 
  std::vector<VBOIndexedTri> fluid_facevelocity_tri_indices;
  std::vector<VBOPosNormalColor> fluid_pressure_vis;
  std::vector<VBOPosNormalColor> fluid_cell_type_vis;

  // Helper class to display an isosurface 
  MarchingCubes *marchingCubes; 
};


void setupCubeVBO(const glm::vec3 pts[8], const glm::vec3 &color, std::vector<VBOPosNormalColor> &faces);

// ========================================================================

#endif
