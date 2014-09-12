#include "glCanvas.h"
#include "marching_cubes.h"
#include "utils.h"

// ============================================================================
// ============================================================================

void MarchingCubes::initializeVBOs() {
  // create a pointer for the vertex & index VBOs
  glGenBuffers(1, &marching_cubes_verts_VBO);
  glGenBuffers(1, &marching_cubes_tri_indices_VBO);
}

void MarchingCubes::setupVBOs() {
  double isosurface = 0.5;
  marching_cubes_verts.clear();
  marching_cubes_tri_indices.clear();

  for (int i = 0; i < nx-1; i++) {
    for (int j = 0; j < ny-1; j++) {
      for (int k = 0; k < nz-1; k++) {
	GridValue v[8];
	double eps = 0;//dx * 0.05;
	v[0] = GridValue(glm::vec3(dx*i    +eps,dy*j    +eps,dz*k    +eps),getNormal(i  ,j  ,k  ),get(i  ,j  ,k  ));
	v[1] = GridValue(glm::vec3(dx*i    +eps,dy*j    +eps,dz*(k+1)-eps),getNormal(i  ,j  ,k+1),get(i  ,j  ,k+1));
	v[2] = GridValue(glm::vec3(dx*i    +eps,dy*(j+1)-eps,dz*k    +eps),getNormal(i  ,j+1,k  ),get(i  ,j+1,k  ));
	v[3] = GridValue(glm::vec3(dx*i    +eps,dy*(j+1)-eps,dz*(k+1)-eps),getNormal(i  ,j+1,k+1),get(i  ,j+1,k+1));
	v[4] = GridValue(glm::vec3(dx*(i+1)-eps,dy*j    +eps,dz*k    +eps),getNormal(i+1,j  ,k  ),get(i+1,j  ,k  ));
	v[5] = GridValue(glm::vec3(dx*(i+1)-eps,dy*j    +eps,dz*(k+1)-eps),getNormal(i+1,j  ,k+1),get(i+1,j  ,k+1));
	v[6] = GridValue(glm::vec3(dx*(i+1)-eps,dy*(j+1)-eps,dz*k    +eps),getNormal(i+1,j+1,k  ),get(i+1,j+1,k  ));
	v[7] = GridValue(glm::vec3(dx*(i+1)-eps,dy*(j+1)-eps,dz*(k+1)-eps),getNormal(i+1,j+1,k+1),get(i+1,j+1,k+1));
	// need to alternate orientation of central tetrahedron to ensure that the diagonals line up
	if ((i+j+k)%2) {
	  PaintTetra(v[0],v[5],v[3],v[6],isosurface);
	  PaintTetra(v[0],v[1],v[3],v[5],isosurface);
	  PaintTetra(v[0],v[4],v[5],v[6],isosurface);
	  PaintTetra(v[0],v[2],v[6],v[3],isosurface);
	  PaintTetra(v[3],v[7],v[6],v[5],isosurface);
	} else {
	  PaintTetra(v[2],v[4],v[1],v[7],isosurface);
	  PaintTetra(v[5],v[4],v[7],v[1],isosurface);
	  PaintTetra(v[2],v[3],v[7],v[1],isosurface);
	  PaintTetra(v[4],v[2],v[1],v[0],isosurface);
	  PaintTetra(v[2],v[7],v[6],v[4],isosurface);
	}
      }
    }
  }

  // copy the data to each VBO
  glBindBuffer(GL_ARRAY_BUFFER,marching_cubes_verts_VBO); 
  glBufferData(GL_ARRAY_BUFFER,sizeof(VBOPosNormalColor)*marching_cubes_verts.size(),&marching_cubes_verts[0],GL_STATIC_DRAW); 
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,marching_cubes_tri_indices_VBO); 
  glBufferData(GL_ELEMENT_ARRAY_BUFFER,sizeof(VBOIndexedTri)*marching_cubes_tri_indices.size(),&marching_cubes_tri_indices[0],GL_STATIC_DRAW);
}


void MarchingCubes::drawVBOs() {
  glUniform1i(GLCanvas::colormodeID, 1);
  glBindBuffer(GL_ARRAY_BUFFER,marching_cubes_verts_VBO); 
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,marching_cubes_tri_indices_VBO); 
  glEnableVertexAttribArray(0);
  glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,3*sizeof(glm::vec3),(void*)0);
  glEnableVertexAttribArray(1);
  glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE,3*sizeof(glm::vec3),(void*)sizeof(glm::vec3) );
  glEnableVertexAttribArray(2);
  glVertexAttribPointer(2, 3, GL_FLOAT,GL_FALSE,3*sizeof(glm::vec3), (void*)(sizeof(glm::vec3)*2));
  glDrawElements(GL_TRIANGLES,
                 marching_cubes_tri_indices.size()*3,
                 GL_UNSIGNED_INT, BUFFER_OFFSET(0));
  glDisableVertexAttribArray(0);
  glDisableVertexAttribArray(1);
  glDisableVertexAttribArray(2);
}


void MarchingCubes::cleanupVBOs() {
  glDeleteBuffers(1, &marching_cubes_verts_VBO);
  glDeleteBuffers(1, &marching_cubes_tri_indices_VBO);
}


// ============================================================================

void MarchingCubes::PaintTetra(const GridValue &a, const GridValue &b, const GridValue &c, const GridValue &d, double isosurface)  {
  // figure out how many of the vertices are on each side of the isosurface
  int count = 0;
  if (a.value > isosurface) count++;
  if (b.value > isosurface) count++;
  if (c.value > isosurface) count++;
  if (d.value > isosurface) count++;
  // handle each case
  if (count == 4) { PaintTetraHelper4(a,b,c,d,isosurface); } 
  else if (count == 3) {
    if (d.value <= isosurface) PaintTetraHelper3(a,b,c,d,isosurface);
    else if (c.value <= isosurface) PaintTetraHelper3(a,d,b,c,isosurface);
    else if (b.value <= isosurface) PaintTetraHelper3(a,c,d,b,isosurface);
    else if (a.value <= isosurface) PaintTetraHelper3(b,d,c,a,isosurface);
    else assert(0);
  } else if (count == 2) {
    if (a.value > isosurface && b.value > isosurface) PaintTetraHelper2(a,b,c,d,isosurface);
    else if (a.value > isosurface && c.value > isosurface) PaintTetraHelper2(a,c,d,b,isosurface);
    else if (a.value > isosurface && d.value > isosurface) PaintTetraHelper2(a,d,b,c,isosurface);
    else if (b.value > isosurface && c.value > isosurface) PaintTetraHelper2(b,c,a,d,isosurface);
    else if (b.value > isosurface && d.value > isosurface) PaintTetraHelper2(b,d,c,a,isosurface);
    else if (c.value > isosurface && d.value > isosurface) PaintTetraHelper2(c,d,a,b,isosurface);
    else assert(0);
  } else if (count == 1) {
    if (a.value > isosurface) PaintTetraHelper1(a,b,c,d,isosurface);
    else if (b.value > isosurface) PaintTetraHelper1(b,c,a,d,isosurface);
    else if (c.value > isosurface) PaintTetraHelper1(c,a,b,d,isosurface);
    else if (d.value > isosurface) PaintTetraHelper1(d,c,b,a,isosurface);
    else assert(0);
  } else {
    // count == 0, do nothing
    assert (count == 0);
  }
}

// ============================================================================

// when all 4 vertices are within the volume
void MarchingCubes::PaintTetraHelper4(const GridValue &a, const GridValue &b, const GridValue &c, const GridValue &d, double isosurface)  {
  assert (a.value > isosurface && b.value > isosurface && c.value > isosurface && d.value > isosurface);
  drawIfBoundary(a.position,b.position,c.position);
  drawIfBoundary(a.position,d.position,b.position);
  drawIfBoundary(b.position,d.position,c.position);
  drawIfBoundary(c.position,d.position,a.position);
}

// when 3 vertices are within the volume
void MarchingCubes::PaintTetraHelper3(const GridValue &a, const GridValue &b, const GridValue &c, const GridValue &d, double isosurface)  {
  assert (a.value > isosurface && b.value > isosurface && c.value > isosurface && d.value <= isosurface);
  float ad_frac = (isosurface - d.value) / (a.value - d.value);
  assert (ad_frac >= 0 && ad_frac <= 1);
  float bd_frac = (isosurface - d.value) / (b.value - d.value);
  assert (bd_frac >= 0 && bd_frac <= 1);
  float cd_frac = (isosurface - d.value) / (c.value - d.value);
  assert (cd_frac >= 0 && cd_frac <= 1);
  glm::vec3 ad = ad_frac*a.position + (1-ad_frac)*d.position;
  glm::vec3 bd = bd_frac*b.position + (1-bd_frac)*d.position;
  glm::vec3 cd = cd_frac*c.position + (1-cd_frac)*d.position;
  glm::vec3 ad_norm = ad_frac*a.normal + (1-ad_frac)*d.normal;
  glm::vec3 bd_norm = bd_frac*b.normal + (1-bd_frac)*d.normal;
  glm::vec3 cd_norm = cd_frac*c.normal + (1-cd_frac)*d.normal;
  drawTriangleWithNormals(ad_norm,ad,cd_norm,cd,bd_norm,bd);
  drawIfBoundary(a.position,ad,bd);
  drawIfBoundary(a.position,bd,b.position);
  drawIfBoundary(b.position,bd,cd);
  drawIfBoundary(b.position,cd,c.position);
  drawIfBoundary(c.position,cd,ad);
  drawIfBoundary(c.position,ad,a.position);
  drawIfBoundary(a.position,b.position,c.position);
}

// when 2 vertices are within the volume
void MarchingCubes::PaintTetraHelper2(const GridValue &a, const GridValue &b, const GridValue &c, const GridValue &d, double isosurface)  {
  assert (a.value > isosurface && b.value > isosurface && c.value <= isosurface && d.value <= isosurface);
  float ac_frac = (isosurface - c.value) / (a.value - c.value);
  assert (ac_frac >= 0 && ac_frac <= 1);
  float ad_frac = (isosurface - d.value) / (a.value - d.value);
  assert (ad_frac >= 0 && ad_frac <= 1);
  float bc_frac = (isosurface - c.value) / (b.value - c.value);
  assert (bc_frac >= 0 && bc_frac <= 1);
  float bd_frac = (isosurface - d.value) / (b.value - d.value);
  assert (bd_frac >= 0 && bd_frac <= 1);
  glm::vec3 ac = ac_frac*a.position + (1-ac_frac)*c.position;
  glm::vec3 ad = ad_frac*a.position + (1-ad_frac)*d.position;
  glm::vec3 bc = bc_frac*b.position + (1-bc_frac)*c.position;
  glm::vec3 bd = bd_frac*b.position + (1-bd_frac)*d.position;
  glm::vec3 ac_norm = ac_frac*a.normal + (1-ac_frac)*c.normal;
  glm::vec3 ad_norm = ad_frac*a.normal + (1-ad_frac)*d.normal;
  glm::vec3 bc_norm = bc_frac*b.normal + (1-bc_frac)*c.normal;
  glm::vec3 bd_norm = bd_frac*b.normal + (1-bd_frac)*d.normal;
  drawTriangleWithNormals(ac_norm,ac,bc_norm,bc,ad_norm,ad);
  drawTriangleWithNormals(bc_norm,bc,bd_norm,bd,ad_norm,ad);
  drawIfBoundary(a.position,ac,ad);
  drawIfBoundary(b.position,bd,bc);
  drawIfBoundary(a.position,ad,bd);
  drawIfBoundary(a.position,bd,b.position);
  drawIfBoundary(a.position,bc,ac);
  drawIfBoundary(a.position,b.position,bc);
}

// when 1 vertex is within the volume
void MarchingCubes::PaintTetraHelper1(const GridValue &a, const GridValue &b, const GridValue &c, const GridValue &d, double isosurface)  {
  assert (a.value > isosurface && b.value <= isosurface && c.value <= isosurface && d.value <= isosurface);
  float ab_frac = (isosurface - b.value) / (a.value - b.value);
  assert (ab_frac >= 0 && ab_frac <= 1);
  float ac_frac = (isosurface - c.value) / (a.value - c.value);
  assert (ac_frac >= 0 && ac_frac <= 1);
  float ad_frac = (isosurface - d.value) / (a.value - d.value);
  assert (ad_frac >= 0 && ad_frac <= 1);
  glm::vec3 ab = ab_frac*a.position + (1-ab_frac)*b.position;
  glm::vec3 ac = ac_frac*a.position + (1-ac_frac)*c.position;
  glm::vec3 ad = ad_frac*a.position + (1-ad_frac)*d.position;
  glm::vec3 ab_norm = ab_frac*a.normal + (1-ab_frac)*b.normal;
  glm::vec3 ac_norm = ac_frac*a.normal + (1-ac_frac)*c.normal;
  glm::vec3 ad_norm = ad_frac*a.normal + (1-ad_frac)*d.normal;
  drawTriangleWithNormals(ab_norm,ab,ad_norm,ad,ac_norm,ac);
  drawIfBoundary(a.position,ab,ac);
  drawIfBoundary(a.position,ac,ad);
  drawIfBoundary(a.position,ad,ab);
}

// ============================================================================

void MarchingCubes::drawTriangleWithNormals(const glm::vec3 &n1, const glm::vec3 &p1, 
					    const glm::vec3 &n2, const glm::vec3 &p2, 
					    const glm::vec3 &n3, const glm::vec3 &p3) {
  int num_marching_cubes_tris = marching_cubes_tri_indices.size();
  glm::vec3 color(0,0,1);
  marching_cubes_verts.push_back(VBOPosNormalColor(p1,n1,color));
  marching_cubes_verts.push_back(VBOPosNormalColor(p2,n2,color));
  marching_cubes_verts.push_back(VBOPosNormalColor(p3,n3,color));
  marching_cubes_tri_indices.push_back(VBOIndexedTri(num_marching_cubes_tris*3+0,
						     num_marching_cubes_tris*3+1,
						     num_marching_cubes_tris*3+2));
}

void MarchingCubes::drawIfBoundary(const glm::vec3 &p1, const glm::vec3 &p2, const glm::vec3 &p3) {
  if ((p1.x <= 0.1*dx && p2.x <= 0.1*dx && p3.x <= 0.1*dx) ||
      (p1.y <= 0.1*dy && p2.y <= 0.1*dy && p3.y <= 0.1*dy) ||
      (p1.z <= 0.1*dz && p2.z <= 0.1*dz && p3.z <= 0.1*dz) ||
      (p1.x >= dx*(nx-1.1) && p2.x >= dx*(nx-1.1) && p3.x >= dx*(nx-1.1)) ||
      (p1.y >= dy*(ny-1.1) && p2.y >= dy*(ny-1.1) && p3.y >= dy*(ny-1.1)) ||
      (p1.z >= dz*(nz-1.1) && p2.z >= dz*(nz-1.1) && p3.z >= dz*(nz-1.1))) {
    glm::vec3 normal = computeNormal(p1,p2,p3);
    drawTriangleWithNormals(normal,p1,normal,p2,normal,p3);
  }
}

// ============================================================================
