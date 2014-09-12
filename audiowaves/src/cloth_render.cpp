#include "glCanvas.h"

#include <fstream>
#include "cloth.h"
#include "argparser.h"
#include "utils.h"


// ================================================================================

void Cloth::initializeVBOs() {
  glGenBuffers(1, &cloth_verts_VBO);
  glGenBuffers(1, &cloth_tri_indices_VBO);

  glGenBuffers(1, &cloth_velocity_verts_VBO);
  glGenBuffers(1, &cloth_velocity_tri_indices_VBO);

  glGenBuffers(1, &cloth_force_verts_VBO);
  glGenBuffers(1, &cloth_force_tri_indices_VBO);
}


void Cloth::AddWireFrameTriangle(const glm::vec3 &apos, const glm::vec3 &bpos, const glm::vec3 &cpos,
                                 const glm::vec3 &anormal, const glm::vec3 &bnormal, const glm::vec3 &cnormal,
                                 const glm::vec3 &abcolor, const glm::vec3 &bccolor, const glm::vec3 &cacolor) {
  glm::vec3 white = glm::vec3(1,1,1);
  glm::vec3 xpos = (apos+bpos+cpos) * (1/3.0f);
  glm::vec3 xnormal = (anormal+bnormal+cnormal);
  xnormal = glm::normalize(xnormal);
  int start = cloth_verts.size();
  cloth_verts.push_back(VBOPosNormalColor(apos,anormal,abcolor));
  cloth_verts.push_back(VBOPosNormalColor(bpos,bnormal,abcolor));
  cloth_verts.push_back(VBOPosNormalColor(xpos,xnormal,white));
  cloth_tri_indices.push_back(VBOIndexedTri(start,start+1,start+2));
  cloth_verts.push_back(VBOPosNormalColor(bpos,bnormal,bccolor));
  cloth_verts.push_back(VBOPosNormalColor(cpos,cnormal,bccolor));
  cloth_verts.push_back(VBOPosNormalColor(xpos,xnormal,white));
  cloth_tri_indices.push_back(VBOIndexedTri(start+3,start+4,start+5));
  cloth_verts.push_back(VBOPosNormalColor(cpos,cnormal,cacolor));
  cloth_verts.push_back(VBOPosNormalColor(apos,anormal,cacolor));
  cloth_verts.push_back(VBOPosNormalColor(xpos,xnormal,white));
  cloth_tri_indices.push_back(VBOIndexedTri(start+6,start+7,start+8));
}


glm::vec3 super_elastic_color(const ClothParticle &a, const ClothParticle &b, double correction) {
  glm::vec3 a_o, b_o, a_c, b_c;
  a_c = a.getPosition();
  b_c = b.getPosition();
  a_o = a.getOriginalPosition();
  b_o = b.getOriginalPosition();
  double length_o,length;
  length = glm::length(a_c-b_c);
  length_o = glm::length(a_o-b_o);

  if (length >= (1+0.99*correction) * length_o){
    // spring is too long, make it cyan
    return glm::vec3(0,1,1);
  } else if (length <= (1-0.99*correction) * length_o) {
    // spring is too short, make it yellow
    return glm::vec3(1,1,0);
  } else {
    return glm::vec3(0,0,0);
  }
}


void Cloth::setupVBOs() {

  HandleGLError("enter setup cloth VBOs");

  cloth_verts.clear();
  cloth_tri_indices.clear();
  cloth_velocity_verts.clear();
  cloth_velocity_tri_indices.clear();
  cloth_force_verts.clear();
  cloth_force_tri_indices.clear();


  // like the last assignment...  to make wireframe edges...
  //
  //   a-----------------------b
  //   |\                     /|
  //   |  \                 /  |
  //   |    \             /    |
  //   |      \         /      |
  //   |        \     /        |
  //   |          \ /          |
  //   |           x           |
  //   |          / \          |
  //   |        /     \        |
  //   |      /         \      |
  //   |    /             \    |
  //   |  /                 \  |
  //   |/                     \|
  //   d-----------------------c
  //


  // mesh surface positions & normals
  for (int i = 0; i < nx-1; i++) {
    for (int j = 0; j < ny-1; j++) {
      
      const ClothParticle &a = getParticle(i,j);
      const ClothParticle &b = getParticle(i,j+1);
      const ClothParticle &c = getParticle(i+1,j+1);
      const ClothParticle &d = getParticle(i+1,j);

      const glm::vec3 &a_pos = a.getPosition();
      const glm::vec3 &b_pos = b.getPosition();
      const glm::vec3 &c_pos = c.getPosition();
      const glm::vec3 &d_pos = d.getPosition();

      glm::vec3 x_pos = (a_pos+b_pos+c_pos+d_pos) * 0.25f;

      glm::vec3 a_normal = computeGouraudNormal(i,j); 
      glm::vec3 b_normal = computeGouraudNormal(i,j+1); 
      glm::vec3 c_normal = computeGouraudNormal(i+1,j+1); 
      glm::vec3 d_normal = computeGouraudNormal(i+1,j); 

      glm::vec3 x_normal = (a_normal+b_normal+c_normal+d_normal);
      x_normal = glm::normalize(x_normal);

      glm::vec3 white = glm::vec3(1,1,1);
      glm::vec3 black = glm::vec3(0,0,0);
      glm::vec3 cyan  = glm::vec3(0,1,1);

      glm::vec3 ab_color = super_elastic_color(a,b,provot_structural_correction);
      glm::vec3 bc_color = super_elastic_color(b,c,provot_structural_correction);
      glm::vec3 cd_color = super_elastic_color(c,d,provot_structural_correction);
      glm::vec3 da_color = super_elastic_color(d,a,provot_structural_correction);

      glm::vec3 ac_color = super_elastic_color(a,c,provot_shear_correction);
      glm::vec3 bd_color = super_elastic_color(b,d,provot_shear_correction);

      AddWireFrameTriangle(a_pos,b_pos,x_pos, 
                           a_normal,b_normal,x_normal,
                           ab_color,bd_color,ac_color);
      AddWireFrameTriangle(b_pos,c_pos,x_pos, 
                           b_normal,c_normal,x_normal,
                           bc_color,ac_color,bd_color);
      AddWireFrameTriangle(c_pos,d_pos,x_pos, 
                           c_normal,d_normal,x_normal,
                           cd_color,bd_color,ac_color);
      AddWireFrameTriangle(d_pos,a_pos,x_pos, 
                           d_normal,a_normal,x_normal,
                           da_color,ac_color,bd_color);

    }
  }


  float thickness = 0.003 * GLCanvas::bbox.maxDim();

  // velocity & force visualization
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      const ClothParticle &p = getParticle(i,j);
      const glm::vec3 &pos = p.getPosition();
      const glm::vec3 &vel = p.getVelocity();

      float dt = args->timestep;

      if (args->velocity) {
        addEdgeGeometry(cloth_velocity_verts,cloth_velocity_tri_indices,
                      pos,pos+dt*100*vel,
                      glm::vec3(1,0,0),glm::vec3(1,0,0),thickness,thickness);
      }


      if (args->force) {


      // *********************************************************************  
      // ASSIGNMENT:
      //
      // Visualize the forces
      //
      // *********************************************************************    


      }
    }
  }

  // copy the data to each VBO
  glBindBuffer(GL_ARRAY_BUFFER,cloth_verts_VBO); 
  glBufferData(GL_ARRAY_BUFFER,sizeof(VBOPosNormalColor)*cloth_verts.size(),&cloth_verts[0],GL_STATIC_DRAW); 
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,cloth_tri_indices_VBO); 
  glBufferData(GL_ELEMENT_ARRAY_BUFFER,sizeof(VBOIndexedTri)*cloth_tri_indices.size(),&cloth_tri_indices[0],GL_STATIC_DRAW);

  glBindBuffer(GL_ARRAY_BUFFER,cloth_velocity_verts_VBO); 
  glBufferData(GL_ARRAY_BUFFER,sizeof(VBOPosNormalColor)*cloth_velocity_verts.size(),&cloth_velocity_verts[0],GL_STATIC_DRAW); 
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,cloth_velocity_tri_indices_VBO); 
  glBufferData(GL_ELEMENT_ARRAY_BUFFER,sizeof(VBOIndexedTri)*cloth_velocity_tri_indices.size(),&cloth_velocity_tri_indices[0],GL_STATIC_DRAW);

  glBindBuffer(GL_ARRAY_BUFFER,cloth_force_verts_VBO); 
  glBufferData(GL_ARRAY_BUFFER,sizeof(VBOPosNormalColor)*cloth_force_verts.size(),&cloth_force_verts[0],GL_STATIC_DRAW); 
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,cloth_force_tri_indices_VBO); 
  glBufferData(GL_ELEMENT_ARRAY_BUFFER,sizeof(VBOIndexedTri)*cloth_force_tri_indices.size(),&cloth_force_tri_indices[0],GL_STATIC_DRAW);

  HandleGLError("leaving setup cloth");
}

void Cloth::drawVBOs() {

  HandleGLError("enter cloth::drawVBOs()");

  // =====================================================================================
  // render the particles
  // =====================================================================================
  if (args->particles) {
    glPointSize(10);
    glBindBuffer(GL_ARRAY_BUFFER,cloth_verts_VBO); 
    glEnableVertexAttribArray(0);
    // There are 8 *extra* points for every grid cell, so set the stride to be 9 vertices....
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(VBOPosNormalColor)*9, 0);
    glDrawArrays(GL_POINTS, 0, cloth_verts.size()/9);
    glDisableVertexAttribArray(0);
  }
  HandleGLError("points done cloth::drawVBOs()");

  // =====================================================================================
  // render the cloth surface
  // =====================================================================================
  if (args->surface) {
    glUniform1i(GLCanvas::colormodeID, 0);
    glBindBuffer(GL_ARRAY_BUFFER,cloth_verts_VBO); 
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,cloth_tri_indices_VBO); 
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,3*sizeof(glm::vec3),(void*)0);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE,3*sizeof(glm::vec3),(void*)sizeof(glm::vec3) );
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 3, GL_FLOAT,GL_FALSE,3*sizeof(glm::vec3), (void*)(sizeof(glm::vec3)*2));
    glDrawElements(GL_TRIANGLES,
                   cloth_tri_indices.size()*3,
                   GL_UNSIGNED_INT, BUFFER_OFFSET(0));
    glDisableVertexAttribArray(0);
    glDisableVertexAttribArray(1);
    glDisableVertexAttribArray(2);
  }
  
  // =====================================================================================
  // render the velocity at each particle
  // =====================================================================================
  if (args->velocity) {
    glUniform1i(GLCanvas::colormodeID, 1);
    glBindBuffer(GL_ARRAY_BUFFER,cloth_velocity_verts_VBO); 
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,cloth_velocity_tri_indices_VBO); 
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,3*sizeof(glm::vec3),(void*)0);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE,3*sizeof(glm::vec3),(void*)sizeof(glm::vec3) );
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 3, GL_FLOAT,GL_FALSE,3*sizeof(glm::vec3), (void*)(sizeof(glm::vec3)*2));
    glDrawElements(GL_TRIANGLES,
                   cloth_velocity_tri_indices.size()*3,
                   GL_UNSIGNED_INT, BUFFER_OFFSET(0));
    glDisableVertexAttribArray(0);
    glDisableVertexAttribArray(1);
    glDisableVertexAttribArray(2);
  }
    
  // =====================================================================================
  // render the forces at each particle
  // =====================================================================================
  if (args->force) {

    
    // *********************************************************************  
    // ASSIGNMENT:
    //
    // Implement this visualization.
    //
    // *********************************************************************    



  }

  HandleGLError("leaving cloth::drawVBOs()");
}

void Cloth::cleanupVBOs() { 
  glDeleteBuffers(1, &cloth_verts_VBO);
  glDeleteBuffers(1, &cloth_tri_indices_VBO);
 
  glDeleteBuffers(1, &cloth_velocity_verts_VBO);
  glDeleteBuffers(1, &cloth_velocity_tri_indices_VBO);

  glDeleteBuffers(1, &cloth_force_verts_VBO);
  glDeleteBuffers(1, &cloth_force_tri_indices_VBO);
}



// ================================================================================
// some helper functions
// ================================================================================


glm::vec3 Cloth::computeGouraudNormal(int i, int j) const {
  assert (i >= 0 && i < nx && j >= 0 && j < ny);

  glm::vec3 pos = getParticle(i,j).getPosition();
  glm::vec3 north = pos;
  glm::vec3 south = pos;
  glm::vec3 east = pos;
  glm::vec3 west = pos;
  
  if (i-1 >= 0) north = getParticle(i-1,j).getPosition();
  if (i+1 < nx) south = getParticle(i+1,j).getPosition();
  if (j-1 >= 0) east = getParticle(i,j-1).getPosition();
  if (j+1 < ny) west = getParticle(i,j+1).getPosition();

  glm::vec3 vns = north - south;
  glm::vec3 vwe = west - east;
  vns = glm::normalize(vns);
  vwe = glm::normalize(vwe);

  // compute normals at each corner and average them
  glm::vec3 normal = glm::cross(vns,vwe);
  normal = glm::normalize(normal);
  return normal;
}

// ================================================================================
