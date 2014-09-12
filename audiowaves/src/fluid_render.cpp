#include "glCanvas.h"

#include <fstream>
#include <algorithm>

#include "fluid.h"
#include "argparser.h"
#include "boundingbox.h"
#include "marching_cubes.h"
#include "utils.h"

// ==============================================================
// ==============================================================

void Fluid::initializeVBOs() {
  glGenBuffers(1, &fluid_particles_VBO);
  glGenBuffers(1, &fluid_velocity_verts_VBO);
  glGenBuffers(1, &fluid_velocity_tri_indices_VBO);
  glGenBuffers(1, &fluid_facevelocity_verts_VBO);
  glGenBuffers(1, &fluid_facevelocity_tri_indices_VBO);
  glGenBuffers(1, &fluid_pressure_vis_VBO);
  glGenBuffers(1, &fluid_cell_type_vis_VBO);
  marchingCubes->initializeVBOs();
}


void Fluid::setupVBOs() {
  HandleGLError("in setup fluid VBOs");

  fluid_particles.clear();
  fluid_velocity_verts.clear();
  fluid_velocity_tri_indices.clear();
  fluid_facevelocity_verts.clear();
  fluid_facevelocity_tri_indices.clear();
  fluid_pressure_vis.clear();
  fluid_cell_type_vis.clear();

  float thickness = 0.002 * GLCanvas::bbox.maxDim();

  // =====================================================================================
  // setup the particles
  // =====================================================================================
  for (int x = 0; x < nx; x++) {
    for (int y = 0; y < ny; y++) {
      for (int z = 0; z < nz; z++) {
	Cell *cell = getCell(x,y,z);
	std::vector<FluidParticle*> &particles = cell->getParticles();
	for (unsigned int iter = 0; iter < particles.size(); iter++) {
	  FluidParticle *p = particles[iter];
	  glm::vec3 v = p->getPosition();
          fluid_particles.push_back(VBOPosNormalColor(v));
	}
      }
    }
  }

  // =====================================================================================
  // visualize the velocity
  // =====================================================================================
  if (args->velocity) {
    if (args->dense_velocity == 0) {
      // one velocity vector per cell, at the centroid
      for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
          for (int k = 0; k < nz; k++) {
            glm::vec3 cell_center((i+0.5)*dx,(j+0.5)*dy,(k+0.5)*dz);
            glm::vec3 direction(get_u_avg(i,j,k),get_v_avg(i,j,k),get_w_avg(i,j,k));
            glm::vec3 pt2 = cell_center+float(100*args->timestep)*direction;
            addEdgeGeometry(fluid_velocity_verts,fluid_velocity_tri_indices,
                            cell_center,pt2,
                            glm::vec3(1,0,0),glm::vec3(1,0,0),thickness,thickness*0.1);
          }
        }
      }
    } else if (args->dense_velocity == 1) {
      double z = nz*dz / 2.0;
      for (double x = 0; x <= (nx+0.01)*dx; x+=0.25*dx) {
        for (double y = 0; y <= (ny+0.01)*dy; y+=0.25*dy) {
          glm::vec3 vel = getInterpolatedVelocity(glm::vec3(x,y,z));
          glm::vec3 pt1(x,y,z);
          glm::vec3 pt2 = pt1 + float(100*args->timestep)*vel;
          addEdgeGeometry(fluid_velocity_verts,fluid_velocity_tri_indices,
                          pt1,pt2,
                          glm::vec3(1,0,0),glm::vec3(1,0,0),thickness,thickness*0.1);
        } 
      }
    } else if (args->dense_velocity == 2) {
      double y = ny*dy / 2.0;
      for (double x = 0; x <= (nx+0.01)*dx; x+=0.25*dx) {
        for (double z = 0; z <= (nz+0.01)*dz; z+=0.25*dz) {
          glm::vec3 vel = getInterpolatedVelocity(glm::vec3(x,y,z));
          glm::vec3 pt1(x,y,z);
          glm::vec3 pt2 = pt1 + float(100*args->timestep)*vel;
          addEdgeGeometry(fluid_velocity_verts,fluid_velocity_tri_indices,
                          pt1,pt2,
                          glm::vec3(1,0,0),glm::vec3(1,0,0),thickness,thickness*0.1);
        }
      } 
    } else if (args->dense_velocity == 3) {
      double x = nx*dx / 2.0;
      for (double y = 0; y <= (ny+0.01)*dy; y+=0.25*dy) {
        for (double z = 0; z <= (nz+0.01)*dz; z+=0.25*dz) {
          glm::vec3 vel = getInterpolatedVelocity(glm::vec3(x,y,z));
          glm::vec3 pt1(x,y,z);
          glm::vec3 pt2 = pt1 + float(100*args->timestep)*vel;
          addEdgeGeometry(fluid_velocity_verts,fluid_velocity_tri_indices,
                          pt1,pt2,
                          glm::vec3(1,0,0),glm::vec3(1,0,0),thickness,thickness*0.1);
        }
      } 
    }
  }


  // =====================================================================================
  // visualize the face velocity
  // render stubby triangles to visualize the u, v, and w velocities between cell faces
  // =====================================================================================
  if (args->face_velocity) {
    for (int i = 0; i < nx; i++) {
      for (int j = 0; j < ny; j++) {
        for (int k = 0; k < nz; k++) {
          double dt = args->timestep;
          double u = get_u_plus(i,j,k)*100*dt;
          double v = get_v_plus(i,j,k)*100*dt;
          double w = get_w_plus(i,j,k)*100*dt;
          double x = i*dx;
          double y = j*dy;
          double z = k*dz;
          addEdgeGeometry(fluid_facevelocity_verts,fluid_facevelocity_tri_indices,
                          glm::vec3(x+dx  ,y+0.5*dy,z+0.5*dz),
                          glm::vec3(x+dx+u,y+0.5*dy,z+0.5*dz),
                          glm::vec3(1,0,0),glm::vec3(1,0,0),thickness,thickness*0.1);
          addEdgeGeometry(fluid_facevelocity_verts,fluid_facevelocity_tri_indices,
                          glm::vec3(x+0.5*dx,y+dy,z+0.5*dz),
                          glm::vec3(x+0.5*dx,y+dy+v,z+0.5*dz),
                          glm::vec3(0,1,0),glm::vec3(0,1,0),thickness,thickness*0.1);
          addEdgeGeometry(fluid_facevelocity_verts,fluid_facevelocity_tri_indices,
                          glm::vec3(x+0.5*dx,y+0.5*dy,z+dz),
                          glm::vec3(x+0.5*dx,y+0.5*dy,z+dz+w),
                          glm::vec3(0,0,1),glm::vec3(0,0,1),thickness,thickness*0.1);
        }
      }
    }
  }

    
  // =====================================================================================
  // visualize the cell pressure
  // =====================================================================================
  if (args->pressure) {
    for (int i = 0; i < nx; i++) {
      for (int j = 0; j < ny; j++) {
        for (int k = 0; k < nz; k++) {
          glm::vec3 pts[8] = { glm::vec3((i+0.1)*dx,(j+0.1)*dy,(k+0.1)*dz),
                               glm::vec3((i+0.1)*dx,(j+0.1)*dy,(k+0.9)*dz),
                               glm::vec3((i+0.1)*dx,(j+0.9)*dy,(k+0.1)*dz),
                               glm::vec3((i+0.1)*dx,(j+0.9)*dy,(k+0.9)*dz),
                               glm::vec3((i+0.9)*dx,(j+0.1)*dy,(k+0.1)*dz),
                               glm::vec3((i+0.9)*dx,(j+0.1)*dy,(k+0.9)*dz),
                               glm::vec3((i+0.9)*dx,(j+0.9)*dy,(k+0.1)*dz),
                               glm::vec3((i+0.9)*dx,(j+0.9)*dy,(k+0.9)*dz) };
          double p = getCell(i,j,k)->getPressure();
          // scale the pressure
          p *= 0.01;
          if (p > 1) p = 1;
          if (p < -1) p = -1;
          assert(p >= -1 && p <= 1);
          glm::vec3 color;
          if (p < 0) {
            // negative pressure is blue
            color = glm::vec3(1+p,1+p,1);
          } else {
            // positive pressure is red
            color = glm::vec3(1,1-p,1-p);
          }
          setupCubeVBO(pts,color,fluid_pressure_vis);
        }
      }
    }
  }


  // =====================================================================================
  // render the MAC cells (FULL, SURFACE, or EMPTY)
  // =====================================================================================
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      for (int k = 0; k < nz; k++) {
	glm::vec3 pts[8] = { glm::vec3((i+0.1)*dx,(j+0.1)*dy,(k+0.1)*dz),
			 glm::vec3((i+0.1)*dx,(j+0.1)*dy,(k+0.9)*dz),
			 glm::vec3((i+0.1)*dx,(j+0.9)*dy,(k+0.1)*dz),
			 glm::vec3((i+0.1)*dx,(j+0.9)*dy,(k+0.9)*dz),
			 glm::vec3((i+0.9)*dx,(j+0.1)*dy,(k+0.1)*dz),
			 glm::vec3((i+0.9)*dx,(j+0.1)*dy,(k+0.9)*dz),
			 glm::vec3((i+0.9)*dx,(j+0.9)*dy,(k+0.1)*dz),
			 glm::vec3((i+0.9)*dx,(j+0.9)*dy,(k+0.9)*dz) };
	Cell *cell = getCell(i,j,k);
	glm::vec3 color;
	if (cell->getStatus() == CELL_FULL) {
	  color = glm::vec3(1,0,0);
	} else if (cell->getStatus() == CELL_SURFACE) {
	  color=glm::vec3(0,0,1);
	} else {
	  continue;
	}
	setupCubeVBO(pts,color,fluid_cell_type_vis);
      }
    }
  }

  // copy the data to each VBO
  glBindBuffer(GL_ARRAY_BUFFER,fluid_particles_VBO); 
  glBufferData(GL_ARRAY_BUFFER,sizeof(VBOPosNormalColor)*fluid_particles.size(),&fluid_particles[0],GL_STATIC_DRAW); 

  glBindBuffer(GL_ARRAY_BUFFER,fluid_velocity_verts_VBO); 
  glBufferData(GL_ARRAY_BUFFER,sizeof(VBOPosNormalColor)*fluid_velocity_verts.size(),&fluid_velocity_verts[0],GL_STATIC_DRAW); 
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,fluid_velocity_tri_indices_VBO); 
  glBufferData(GL_ELEMENT_ARRAY_BUFFER,sizeof(VBOIndexedTri)*fluid_velocity_tri_indices.size(),&fluid_velocity_tri_indices[0],GL_STATIC_DRAW);

  glBindBuffer(GL_ARRAY_BUFFER,fluid_facevelocity_verts_VBO); 
  glBufferData(GL_ARRAY_BUFFER,sizeof(VBOPosNormalColor)*fluid_facevelocity_verts.size(),&fluid_facevelocity_verts[0],GL_STATIC_DRAW); 
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,fluid_facevelocity_tri_indices_VBO); 
  glBufferData(GL_ELEMENT_ARRAY_BUFFER,sizeof(VBOIndexedTri)*fluid_facevelocity_tri_indices.size(),&fluid_facevelocity_tri_indices[0],GL_STATIC_DRAW);

  glBindBuffer(GL_ARRAY_BUFFER,fluid_pressure_vis_VBO); 
  glBufferData(GL_ARRAY_BUFFER,sizeof(VBOPosNormalColor)*fluid_pressure_vis.size(),&fluid_pressure_vis[0],GL_STATIC_DRAW); 
  glBindBuffer(GL_ARRAY_BUFFER,fluid_cell_type_vis_VBO); 
  glBufferData(GL_ARRAY_BUFFER,sizeof(VBOPosNormalColor)*fluid_cell_type_vis.size(),&fluid_cell_type_vis[0],GL_STATIC_DRAW); 

  // =====================================================================================
  // setup a marching cubes representation of the surface
  // =====================================================================================
  for (int i = 0; i <= nx; i++) {
    for (int j = 0; j <= ny; j++) {
      for (int k = 0; k <= nz; k++) {
	marchingCubes->set(i,j,k,interpolateIsovalue(glm::vec3((i-0.5),(j-0.5),(k-0.5))));
      } 
    }
  }
  marchingCubes->setupVBOs();

  HandleGLError("leaving setup fluid");
}


void Fluid::drawVBOs() {

  HandleGLError("enter fluid drawVBOs");
  // =====================================================================================
  // render the particles
  // =====================================================================================
  if (args->particles) {
    glPointSize(5);
    glBindBuffer(GL_ARRAY_BUFFER, fluid_particles_VBO);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,3*sizeof(glm::vec3),(void*)0);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE,3*sizeof(glm::vec3),(void*)sizeof(glm::vec3) );
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 3, GL_FLOAT,GL_FALSE,3*sizeof(glm::vec3), (void*)(sizeof(glm::vec3)*2));
    glDrawArrays(GL_POINTS, 0, fluid_particles.size());
    glDisableVertexAttribArray(0);
    glDisableVertexAttribArray(1);
    glDisableVertexAttribArray(2);
  }

  // =====================================================================================
  // visualize the average cell velocity
  // =====================================================================================
  if (args->velocity) {
    glUniform1i(GLCanvas::colormodeID, 1);
    glBindBuffer(GL_ARRAY_BUFFER,fluid_velocity_verts_VBO); 
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,fluid_velocity_tri_indices_VBO); 
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,3*sizeof(glm::vec3),(void*)0);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE,3*sizeof(glm::vec3),(void*)sizeof(glm::vec3) );
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 3, GL_FLOAT,GL_FALSE,3*sizeof(glm::vec3), (void*)(sizeof(glm::vec3)*2));
    glDrawElements(GL_TRIANGLES,
                   fluid_velocity_tri_indices.size()*3,
                   GL_UNSIGNED_INT, BUFFER_OFFSET(0));
    glDisableVertexAttribArray(0);
    glDisableVertexAttribArray(1);
    glDisableVertexAttribArray(2);
  }

  // =====================================================================================
  // visualize the face velocity
  // =====================================================================================
  if (args->face_velocity) {// && fluid_facevelocity_vis.size() > 0) {
    glUniform1i(GLCanvas::colormodeID, 1);
    glBindBuffer(GL_ARRAY_BUFFER,fluid_facevelocity_verts_VBO); 
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,fluid_facevelocity_tri_indices_VBO); 
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,3*sizeof(glm::vec3),(void*)0);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE,3*sizeof(glm::vec3),(void*)sizeof(glm::vec3) );
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 3, GL_FLOAT,GL_FALSE,3*sizeof(glm::vec3), (void*)(sizeof(glm::vec3)*2));
    glDrawElements(GL_TRIANGLES,
                   fluid_facevelocity_tri_indices.size()*3,
                   GL_UNSIGNED_INT, BUFFER_OFFSET(0));
    glDisableVertexAttribArray(0);
    glDisableVertexAttribArray(1);
    glDisableVertexAttribArray(2);
  }

  // =====================================================================================
  // visualize the cell pressure
  // =====================================================================================
  if (args->pressure) {
    glUniform1i(GLCanvas::colormodeID, 1);
    glBindBuffer(GL_ARRAY_BUFFER, fluid_pressure_vis_VBO);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,3*sizeof(glm::vec3),(void*)0);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE,3*sizeof(glm::vec3),(void*)sizeof(glm::vec3) );
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 3, GL_FLOAT,GL_FALSE,3*sizeof(glm::vec3), (void*)(sizeof(glm::vec3)*2));
    glDrawArrays(GL_TRIANGLES, 0, fluid_pressure_vis.size());
    glDisableVertexAttribArray(0);
    glDisableVertexAttribArray(1);
    glDisableVertexAttribArray(2);
  }

  // =====================================================================================
  // render the MAC cells (FULL, SURFACE, or EMPTY)
  // =====================================================================================
  if (args->cubes) {
    glUniform1i(GLCanvas::colormodeID, 1);
    glBindBuffer(GL_ARRAY_BUFFER, fluid_cell_type_vis_VBO);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,3*sizeof(glm::vec3),(void*)0);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE,3*sizeof(glm::vec3),(void*)sizeof(glm::vec3) );
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 3, GL_FLOAT,GL_FALSE,3*sizeof(glm::vec3), (void*)(sizeof(glm::vec3)*2));
    glDrawArrays(GL_TRIANGLES, 0, fluid_cell_type_vis.size());
    glDisableVertexAttribArray(0);
    glDisableVertexAttribArray(1);
    glDisableVertexAttribArray(2);
  }

  // =====================================================================================
  // render a marching cubes representation of the surface
  //    note: make sure you set the Fluid::getIsovalue() to what you want
  // =====================================================================================
  if (args->surface) {
    marchingCubes->drawVBOs();
  } 
}


void Fluid::cleanupVBOs() { 
  glDeleteBuffers(1, &fluid_particles_VBO);
  glDeleteBuffers(1, &fluid_velocity_verts_VBO);
  glDeleteBuffers(1, &fluid_velocity_tri_indices_VBO);
  glDeleteBuffers(1, &fluid_facevelocity_verts_VBO);
  glDeleteBuffers(1, &fluid_facevelocity_tri_indices_VBO);
  glDeleteBuffers(1, &fluid_pressure_vis_VBO);
  glDeleteBuffers(1, &fluid_cell_type_vis_VBO);
  marchingCubes->cleanupVBOs();
}

// ==============================================================

double Fluid::getIsovalue(int i, int j, int k) const {
  i = my_max(0,(my_min(i,nx-1)));
  j = my_max(0,(my_min(j,ny-1)));
  k = my_max(0,(my_min(k,nz-1)));
  Cell *c = getCell(i,j,k);
  if (c->getStatus() == CELL_EMPTY) return 0;
  // note: this is technically not a correct thing to do
  //       the number of particles is not an indication of it's "fullness"
  if (c->getStatus() == CELL_SURFACE) return 0.5 + c->numParticles()/double(density);
  if (c->getStatus() == CELL_FULL) return 2;
  assert(0);
  return 0;
}

// ==============================================================

double Fluid::interpolateIsovalue(const glm::vec3 &v) const {

  double x = v.x;
  double y = v.y;
  double z = v.z;

  // get the values at the corners
  double a = getIsovalue(int(floor(x)),int(floor(y)),int(floor(z)));
  double b = getIsovalue(int(floor(x)),int(floor(y)),int( ceil(z)));
  double c = getIsovalue(int(floor(x)),int( ceil(y)),int(floor(z)));
  double d = getIsovalue(int(floor(x)),int( ceil(y)),int( ceil(z)));
  double e = getIsovalue(int( ceil(x)),int(floor(y)),int(floor(z)));
  double f = getIsovalue(int( ceil(x)),int(floor(y)),int( ceil(z)));
  double g = getIsovalue(int( ceil(x)),int( ceil(y)),int(floor(z)));
  double h = getIsovalue(int( ceil(x)),int( ceil(y)),int( ceil(z)));

  double x_frac = x - (floor(x));
  double y_frac = y - (floor(y));
  double z_frac = z - (floor(z));

  assert (x_frac >= 0 && x_frac <= 1);
  assert (y_frac >= 0 && y_frac <= 1);
  assert (z_frac >= 0 && z_frac <= 1);
  
  double answer = triInterpolate(x_frac,y_frac,z_frac,a,b,c,d,e,f,g,h);
  return answer;
}


// ==============================================================

void setupCubeVBO(const glm::vec3 pts[8], const glm::vec3 &color, std::vector<VBOPosNormalColor> &faces) {
  
  faces.push_back(VBOPosNormalColor(pts[0],glm::vec3(-1,0,0),color));
  faces.push_back(VBOPosNormalColor(pts[1],glm::vec3(-1,0,0),color));
  faces.push_back(VBOPosNormalColor(pts[3],glm::vec3(-1,0,0),color));
  faces.push_back(VBOPosNormalColor(pts[0],glm::vec3(-1,0,0),color));
  faces.push_back(VBOPosNormalColor(pts[3],glm::vec3(-1,0,0),color));
  faces.push_back(VBOPosNormalColor(pts[2],glm::vec3(-1,0,0),color));
  
  faces.push_back(VBOPosNormalColor(pts[4],glm::vec3(1,0,0),color));
  faces.push_back(VBOPosNormalColor(pts[6],glm::vec3(1,0,0),color));
  faces.push_back(VBOPosNormalColor(pts[7],glm::vec3(1,0,0),color));
  faces.push_back(VBOPosNormalColor(pts[4],glm::vec3(1,0,0),color));
  faces.push_back(VBOPosNormalColor(pts[7],glm::vec3(1,0,0),color));
  faces.push_back(VBOPosNormalColor(pts[5],glm::vec3(1,0,0),color));
  
  faces.push_back(VBOPosNormalColor(pts[0],glm::vec3(0,0,-1),color));
  faces.push_back(VBOPosNormalColor(pts[2],glm::vec3(0,0,-1),color));
  faces.push_back(VBOPosNormalColor(pts[6],glm::vec3(0,0,-1),color));
  faces.push_back(VBOPosNormalColor(pts[0],glm::vec3(0,0,-1),color));
  faces.push_back(VBOPosNormalColor(pts[6],glm::vec3(0,0,-1),color));
  faces.push_back(VBOPosNormalColor(pts[4],glm::vec3(0,0,-1),color));
  
  faces.push_back(VBOPosNormalColor(pts[1],glm::vec3(0,0,1),color));
  faces.push_back(VBOPosNormalColor(pts[5],glm::vec3(0,0,1),color));
  faces.push_back(VBOPosNormalColor(pts[7],glm::vec3(0,0,1),color));
  faces.push_back(VBOPosNormalColor(pts[1],glm::vec3(0,0,1),color));
  faces.push_back(VBOPosNormalColor(pts[7],glm::vec3(0,0,1),color));
  faces.push_back(VBOPosNormalColor(pts[3],glm::vec3(0,0,1),color));
  
  faces.push_back(VBOPosNormalColor(pts[0],glm::vec3(0,-1,0),color));
  faces.push_back(VBOPosNormalColor(pts[4],glm::vec3(0,-1,0),color));
  faces.push_back(VBOPosNormalColor(pts[5],glm::vec3(0,-1,0),color));
  faces.push_back(VBOPosNormalColor(pts[0],glm::vec3(0,-1,0),color));
  faces.push_back(VBOPosNormalColor(pts[5],glm::vec3(0,-1,0),color));
  faces.push_back(VBOPosNormalColor(pts[1],glm::vec3(0,-1,0),color));
  
  faces.push_back(VBOPosNormalColor(pts[2],glm::vec3(0,1,0),color));
  faces.push_back(VBOPosNormalColor(pts[3],glm::vec3(0,1,0),color));
  faces.push_back(VBOPosNormalColor(pts[7],glm::vec3(0,1,0),color));
  faces.push_back(VBOPosNormalColor(pts[2],glm::vec3(0,1,0),color));
  faces.push_back(VBOPosNormalColor(pts[7],glm::vec3(0,1,0),color));
  faces.push_back(VBOPosNormalColor(pts[6],glm::vec3(0,1,0),color));
}

// ==============================================================
