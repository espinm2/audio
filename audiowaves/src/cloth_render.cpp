#include "glCanvas.h"

#include <fstream>
#include <iostream>
#include "cloth.h"
#include "argparser.h"
#include "vectors.h"
#include "utils.h"


// ================================================================================

void Cloth::initializeVBOs() {
  glGenBuffers(1, &cloth_verts_VBO);
  glGenBuffers(1, &cloth_quad_indices_VBO);
  glGenBuffers(1, &cloth_happy_edge_indices_VBO);
  glGenBuffers(1, &cloth_unhappy_edge_indices_VBO);
  glGenBuffers(1, &cloth_velocity_visualization_VBO);
  glGenBuffers(1, &cloth_force_visualization_VBO);
}


void Cloth::setupVBOs() {

  HandleGLError("in setup cloth VBOs");
  cloth_verts.clear();
  cloth_quad_indices.clear();
  cloth_happy_edge_indices.clear();
  cloth_unhappy_edge_indices.clear();
  cloth_velocity_visualization.clear();
  cloth_force_visualization.clear();

  // mesh surface positions & normals
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      const ClothParticle &p = getParticle(i,j);
      const Vec3f &pos = p.getPosition();
      Vec3f normal = computeGouraudNormal(i,j); 
      Vec3f color = Vec3f(0,0,0);
      if (p.isFixed()) 
	color = Vec3f(0,1,0);
      cloth_verts.push_back(VBOPosNormalColor(pos,normal,color));
    }
  }
  // mesh surface
  for (int i = 0; i < nx-1; i++) {
    for (int j = 0; j < ny-1; j++) {
      cloth_quad_indices.push_back(VBOIndexedQuad(i*(ny)+j,i*(ny)+j+1,(i+1)*(ny)+j+1,(i+1)*(ny)+j));
    }
  }

  // spring over-/under-stretch visualization
  for (int i = 0; i < nx-1; i++) {
    for (int j = 0; j < ny-1; j++) {
      AddVBOEdge(i  ,j  ,i+1,j  ,provot_structural_correction);
      AddVBOEdge(i  ,j  ,i  ,j+1,provot_structural_correction);
      AddVBOEdge(i  ,j  ,i+1,j+1,provot_shear_correction);
      AddVBOEdge(i  ,j+1,i+1,j  ,provot_shear_correction);
    }
  }      
  for (int i = 0; i < nx-1; i++) {
    AddVBOEdge(i  ,ny-1,i+1,ny-1,provot_structural_correction);
  }
  for (int j = 0; j < ny-1; j++) {
    AddVBOEdge(nx-1,j  ,nx-1,j+1,provot_structural_correction);
  }

  // velocity & force visualization
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      const ClothParticle &p = getParticle(i,j);
      const Vec3f &pos = p.getPosition();
      const Vec3f &vel = p.getVelocity();
      const Vec3f &frc = p.getForce();

      float dt = args->timestep;

      // Velocity is visualized
      cloth_velocity_visualization.push_back(VBOPosColor(pos,Vec3f(1,0,0)));
      cloth_velocity_visualization.push_back(VBOPosColor(pos+dt*100*vel,Vec3f(1,1,1)));

      // Force is visualized
      cloth_force_visualization.push_back(VBOPosColor(pos,Vec3f(0,0,1)));
      cloth_force_visualization.push_back(VBOPosColor(pos+dt*1000*frc,Vec3f(1,1,1)));
    }
  }

  // cleanup old buffer data (if any)
  cleanupVBOs();
  // copy the data to each VBO
  glBindBuffer(GL_ARRAY_BUFFER,cloth_verts_VBO); 
  glBufferData(GL_ARRAY_BUFFER,sizeof(VBOPosNormalColor)*cloth_verts.size(),&cloth_verts[0],GL_STATIC_DRAW); 
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,cloth_quad_indices_VBO); 
  glBufferData(GL_ELEMENT_ARRAY_BUFFER,sizeof(VBOIndexedQuad)*cloth_quad_indices.size(),&cloth_quad_indices[0],GL_STATIC_DRAW);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,cloth_happy_edge_indices_VBO); 
  glBufferData(GL_ELEMENT_ARRAY_BUFFER,sizeof(VBOIndexedEdge)*cloth_happy_edge_indices.size(),&cloth_happy_edge_indices[0],GL_STATIC_DRAW);
  if (cloth_unhappy_edge_indices.size() > 0) {
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,cloth_unhappy_edge_indices_VBO); 
    glBufferData(GL_ELEMENT_ARRAY_BUFFER,sizeof(VBOIndexedEdge)*cloth_unhappy_edge_indices.size(),&cloth_unhappy_edge_indices[0],GL_STATIC_DRAW);
  }
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,cloth_velocity_visualization_VBO); 
  glBufferData(GL_ELEMENT_ARRAY_BUFFER,sizeof(VBOPosColor)*cloth_velocity_visualization.size(),&cloth_velocity_visualization[0],GL_STATIC_DRAW);

  // *********************************************************************  
  // ASSIGNMENT: uncomment/edit as needed for force visualization
  //
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,cloth_force_visualization_VBO); 
  glBufferData(GL_ELEMENT_ARRAY_BUFFER,sizeof(VBOPosColor)*cloth_force_visualization.size(),&cloth_force_visualization[0],GL_STATIC_DRAW);
  //
  // *********************************************************************    

  HandleGLError("leaving setup cloth");
}

void Cloth::drawVBOs() { 

  // =====================================================================================
  // render the particles
  // =====================================================================================
  if (args->particles) {
    glPointSize(3);
    glColor3f(1,0,0);
    glBindBuffer(GL_ARRAY_BUFFER, cloth_verts_VBO);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(VBOPosNormalColor), 0);
    glEnableClientState(GL_VERTEX_ARRAY);
    glVertexPointer(3, GL_FLOAT,sizeof(VBOPosNormalColor), 0);
    glEnableClientState(GL_COLOR_ARRAY);
    glColorPointer(3, GL_FLOAT, sizeof(VBOPosNormalColor),BUFFER_OFFSET(24));
    glDrawArrays(GL_POINTS, 0, nx*ny);
    glDisableClientState(GL_VERTEX_ARRAY);
    glDisableClientState(GL_COLOR_ARRAY);
    glDisableVertexAttribArray(0);
  }


  // =====================================================================================
  // render the cloth surface
  // =====================================================================================
  if (args->surface) {
    glEnable(GL_LIGHTING);
    glEnable(GL_POLYGON_OFFSET_FILL); 
    glPolygonOffset(1.1,4.0); 
    glColor3f(1,1,1);
    glBindBuffer(GL_ARRAY_BUFFER, cloth_verts_VBO);
    glEnableClientState(GL_VERTEX_ARRAY);
    glVertexPointer(3, GL_FLOAT, sizeof(VBOPosNormalColor), BUFFER_OFFSET(0));
    glEnableClientState(GL_NORMAL_ARRAY);
    glNormalPointer(GL_FLOAT, sizeof(VBOPosNormalColor), BUFFER_OFFSET(12));
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, cloth_quad_indices_VBO);
    glDrawElements(GL_QUADS,(nx-1)*(ny-1)*4,GL_UNSIGNED_INT, BUFFER_OFFSET(0));
    glDisableClientState(GL_NORMAL_ARRAY);
    glDisableClientState(GL_VERTEX_ARRAY);
    glDisable(GL_POLYGON_OFFSET_FILL); 
    glDisable(GL_LIGHTING);
  }


  // =====================================================================================
  // visualize the structural and shear springs
  // =====================================================================================
  if (args->wireframe) {
    glLineWidth(1);
    glBindBuffer(GL_ARRAY_BUFFER, cloth_verts_VBO);
    glEnableClientState(GL_VERTEX_ARRAY);
    glVertexPointer(3, GL_FLOAT, sizeof(VBOPosNormalColor), BUFFER_OFFSET(0));
    glLineWidth(2);
    // draw all the "happy" edges
    glColor3f(0,0,0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, cloth_happy_edge_indices_VBO);
    glDrawElements(GL_LINES, cloth_happy_edge_indices.size()*2, GL_UNSIGNED_INT, 0);
    // draw all the "unhappy" (over- or under-stretched) edges
    if (cloth_unhappy_edge_indices.size() > 0) {
      glLineWidth(3);
      glColor3f(0,1,1);
      glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, cloth_unhappy_edge_indices_VBO);
      glDrawElements(GL_LINES, cloth_unhappy_edge_indices.size()*2, GL_UNSIGNED_INT, 0);
    }
    glDisableClientState(GL_VERTEX_ARRAY);
  }

  // =====================================================================================
  // render the velocity at each particle
  // =====================================================================================
  if (args->velocity) {
    glLineWidth(2);
    glBindBuffer(GL_ARRAY_BUFFER, cloth_velocity_visualization_VBO);
    glEnableClientState(GL_VERTEX_ARRAY);
    glVertexPointer(3, GL_FLOAT,sizeof(VBOPosColor), 0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(VBOPosColor), 0);
    glEnableClientState(GL_COLOR_ARRAY);
    glColorPointer(3, GL_FLOAT, sizeof(VBOPosColor),BUFFER_OFFSET(12));
    glDrawArrays(GL_LINES, 0, nx*ny*2);
    glDisableClientState(GL_VERTEX_ARRAY);
    glDisableClientState(GL_COLOR_ARRAY);
    glDisableVertexAttribArray(0);
  }

  // =====================================================================================
  // render the forces at each particle
  // =====================================================================================
  if (args->force) {

    glLineWidth(2);
    glBindBuffer(GL_ARRAY_BUFFER, cloth_force_visualization_VBO);
    glEnableClientState(GL_VERTEX_ARRAY);
    glVertexPointer(3, GL_FLOAT,sizeof(VBOPosColor), 0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(VBOPosColor), 0);
    glEnableClientState(GL_COLOR_ARRAY);
    glColorPointer(3, GL_FLOAT, sizeof(VBOPosColor),BUFFER_OFFSET(12));
    glDrawArrays(GL_LINES, 0, nx*ny*2);
    glDisableClientState(GL_VERTEX_ARRAY);
    glDisableClientState(GL_COLOR_ARRAY);
    glDisableVertexAttribArray(0);

  }
}

void Cloth::cleanupVBOs() { 
  glDeleteBuffers(1, &cloth_verts_VBO);
  glDeleteBuffers(1, &cloth_quad_indices_VBO);
  glDeleteBuffers(1, &cloth_happy_edge_indices_VBO);
  glDeleteBuffers(1, &cloth_unhappy_edge_indices_VBO);
  glDeleteBuffers(1, &cloth_velocity_visualization_VBO);
}



// ================================================================================
// some helper functions
// ================================================================================

void Cloth::AddVBOEdge(int i1, int j1, int i2, int j2, double correction) {
  Vec3f a_o, b_o, a, b;
  a = getParticle(i1,j1).getPosition();
  b = getParticle(i2,j2).getPosition();
  a_o = getParticle(i1,j1).getOriginalPosition();
  b_o = getParticle(i2,j2).getOriginalPosition();
  double length_o,length;
  length = (a-b).Length();
  length_o = (a_o-b_o).Length();
  if (length >= (1+0.99*correction) * length_o ||
      length <= (1-0.99*correction) * length_o) {
    cloth_unhappy_edge_indices.push_back(VBOIndexedEdge(i1*(ny)+j1,i2*(ny)+j2));
  } else {
    cloth_happy_edge_indices.push_back(VBOIndexedEdge(i1*(ny)+j1,i2*(ny)+j2));
  }  
}



Vec3f Cloth::computeGouraudNormal(int i, int j) const {
  assert (i >= 0 && i < nx && j >= 0 && j < ny);

  Vec3f pos = getParticle(i,j).getPosition();
  Vec3f north = pos;
  Vec3f south = pos;
  Vec3f east = pos;
  Vec3f west = pos;
  
  if (i-1 >= 0) north = getParticle(i-1,j).getPosition();
  if (i+1 < nx) south = getParticle(i+1,j).getPosition();
  if (j-1 >= 0) east = getParticle(i,j-1).getPosition();
  if (j+1 < ny) west = getParticle(i,j+1).getPosition();

  Vec3f vns = north - south;
  Vec3f vwe = west - east;
  vns.Normalize();
  vwe.Normalize();

  // compute normals at each corner and average them
  Vec3f normal;
  Vec3f::Cross3(normal,vns,vwe);
  normal.Normalize();
  return normal;
}

// ================================================================================
