#ifndef _MARCHING_CUBES_
#define _MARCHING_CUBES_

#include <vector>
#include "vbo_structs.h"

// ==================================================================================
// A helper class for the Marching Cubes algorithm

class GridValue {
public:
  GridValue() {}
  GridValue(const glm::vec3 &p, const glm::vec3 &n, double v) : position(p), normal(n), value(v) { } 
  glm::vec3 position;
  glm::vec3 normal;
  double value;
};


// ==================================================================================
// The marching cubes algorithm is used to render an isosurface
// of a signed distance field.

class MarchingCubes {

public:
  MarchingCubes(int _nx, int _ny, int _nz, double _dx, double _dy, double _dz) : 
    nx(_nx), ny(_ny), nz(_nz), dx(_dx), dy(_dy), dz(_dz) { values = new double[nx*ny*nz]; }
  ~MarchingCubes() { delete [] values; }

  // position
  double get(int x, int y, int z) const {
    assert (x >= 0 && x < nx);
    assert (y >= 0 && y < ny);
    assert (z >= 0 && z < nz);
    return values[z*nx*ny + y*nx + x];
  }
  // normal
  glm::vec3 getNormal(int x, int y, int z) const {
    assert (x >= 0 && x < nx);
    assert (y >= 0 && y < ny);
    assert (z >= 0 && z < nz);
    double dx;
    double dy;
    double dz;
    if (x == 0)         dx = get(x+1,y,z)-get(x  ,y,z);
    else if (x == nx-1) dx = get(x  ,y,z)-get(x-1,y,z);
    else                dx = get(x+1,y,z)-get(x-1,y,z);
    if (y == 0)         dy = get(x,y+1,z)-get(x,y  ,z);
    else if (y == ny-1) dy = get(x,y  ,z)-get(x,y-1,z);
    else                dy = get(x,y+1,z)-get(x,y-1,z);
    if (z == 0)         dz = get(x,y,z+1)-get(x,y,z  );
    else if (z == nz-1) dz = get(x,y,z  )-get(x,y,z-1);
    else                dz = get(x,y,z+1)-get(x,y,z-1);
    glm::vec3 norm(-dx,-dy,-dz);
    norm = glm::normalize(norm);
    return norm;
  }
  void set(int x, int y, int z, double v) const {
    assert (x >= 0 && x < nx);
    assert (y >= 0 && y < ny);
    assert (z >= 0 && z < nz);
    values[z*nx*ny + y*nx + x] = v;
  }

  // =============
  // THE DRAW CODE
  void initializeVBOs(); 
  void setupVBOs(); 
  void drawVBOs();
  void cleanupVBOs();

private:

  // private helper functions
  void PaintTetra(const GridValue &a, const GridValue &b, const GridValue &c, const GridValue &d, double isosurface); 
  void PaintTetraHelper4(const GridValue &a, const GridValue &b, const GridValue &c, const GridValue &d, double isosurface); 
  void PaintTetraHelper3(const GridValue &a, const GridValue &b, const GridValue &c, const GridValue &d, double isosurface);
  void PaintTetraHelper2(const GridValue &a, const GridValue &b, const GridValue &c, const GridValue &d, double isosurface); 
  void PaintTetraHelper1(const GridValue &a, const GridValue &b, const GridValue &c, const GridValue &d, double isosurface); 
  void drawIfBoundary(const glm::vec3 &p1, const glm::vec3 &p2, const glm::vec3 &p3);
  void drawTriangleWithNormals(const glm::vec3 &n1, const glm::vec3 &p1, 
			       const glm::vec3 &n2, const glm::vec3 &p2, 
			       const glm::vec3 &n3, const glm::vec3 &p3);

  // ==============
  // REPRESENTATION
  int nx, ny, nz;
  double dx, dy, dz;
  double* values;
  
  GLuint marching_cubes_verts_VBO;
  GLuint marching_cubes_tri_indices_VBO;

  std::vector<VBOPosNormalColor> marching_cubes_verts;
  std::vector<VBOIndexedTri> marching_cubes_tri_indices;
};

// ==================================================================================

#endif

