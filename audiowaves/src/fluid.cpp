#include "glCanvas.h"

#include <fstream>
#include <stdio.h>
#include <iomanip>
#include <algorithm>
#include <math.h>

#include "fluid.h"
#include "argparser.h"
#include "boundingbox.h"
#include "vectors.h"
#include "matrix.h"
#include "marching_cubes.h"
#include "utils.h"

#define BETA_0 1.7
#define EPSILON 0.0001

// ==============================================================
// ==============================================================
// CONSTRUCTOR
// ==============================================================
// ==============================================================

Fluid::Fluid(ArgParser *_args) {
  args = _args;
  Load();
  marchingCubes = new MarchingCubes(nx+1,ny+1,nz+1,dx,dy,dz);
  SetEmptySurfaceFull();
  initializeVBOs();
  setupVBOs();
}

Fluid::~Fluid() { 
  delete [] cells; 
  delete marchingCubes; 
  cleanupVBOs(); 
}

// ==============================================================

void Fluid::Load() {    

  // open the file
  assert (args->fluid_file != "");
  std::ifstream istr(args->fluid_file.c_str());
  assert (istr != NULL);
  std::string token, token2, token3, mode;


  istr >> mode;
  if(mode == "Foster"){
    Foster = true;
    istr >> token;
  }else{
    token = mode;
  }

  // load in the grid size & dimensions
  istr >> nx >> ny >> nz;  assert (token=="grid");
  assert (nx > 0 && ny > 0 && nz > 0);
  istr >> token >> dx >> dy >> dz; assert (token=="cell_dimensions");
  cells = new Cell[(nx+2)*(ny+2)*(nz+2)];

  // simulation parameters
  istr >> token >> token2;  assert (token=="flow");
  if (token2 == "compressible") compressible = true;
  else { assert (token2 == "incompressible"); compressible = false; }
  istr >> token >> token2;  assert (token=="xy_boundary");
  if (token2 == "free_slip") xy_free_slip = true;
  else { assert  (token2 == "no_slip"); xy_free_slip = false; }
  istr >> token >> token2;  assert (token=="yz_boundary");
  if (token2 == "free_slip") yz_free_slip = true;
  else { assert  (token2 == "no_slip"); yz_free_slip = false; }
  istr >> token >> token2;  assert (token=="zx_boundary");
  if (token2 == "free_slip") zx_free_slip = true;
  else { assert  (token2 == "no_slip"); zx_free_slip = false; }
  istr >> token >> viscosity;  assert (token=="viscosity");
  double gravity;
  istr >> token >> gravity;  assert (token=="gravity");
  args->gravity = Vec3f(0,-9.8,0) * gravity;
  
  // initialize marker particles 
  istr >> token >> token2 >> token3;  assert (token=="initial_particles");
  istr >> token >> density;  assert (token=="density");
  GenerateParticles(token2,token3);

  // initialize velocities
  istr >> token >> token2;  assert (token=="initial_velocity");
  if (token2 == "zero") {
    // default is zero
  } else {
    assert (token2 == "random");
    int i,j,k;
    double max_dim = my_max(dx,my_max(dy,dz));
    for (i = -1; i <= nx; i++) {
      for (j = -1; j <= ny; j++) {
        for (k = -1; k <= nz; k++) {
          getCell(i,j,k)->set_u_plus((2*args->mtrand.rand()-1)*max_dim);
    getCell(i,j,k)->set_v_plus((2*args->mtrand.rand()-1)*max_dim);
    getCell(i,j,k)->set_w_plus((2*args->mtrand.rand()-1)*max_dim);
        }
      }
    }
  }
  // read in custom velocities
  while(istr >> token) {
    int i,j,k;
    double velocity;
    assert (token == "u" || token == "v" || token == "w");
    istr >> i >> j >> k >> velocity;
    assert(i >= 0 && i < nx);
    assert(j >= 0 && j < ny);
    assert(k >= 0 && k < nz);
    if      (token == "u") getCell(i,j,k)->set_u_plus(velocity);
    else if (token == "v") getCell(i,j,k)->set_v_plus(velocity);
    else if (token == "w") getCell(i,j,k)->set_w_plus(velocity);
    else assert(0);
  }
  SetBoundaryVelocities();
}

// ==============================================================

bool Fluid::inShape(Vec3f &pos, const std::string &shape) {
  // return true if this point is inside the "shape"
  // defined procedurally (using an implicit surface)
  if (shape == "everywhere") {
    return true;
  } else if (shape == "left") {
    // a blob of particles on the lower left (for the dam)
    return (pos.x() < 0.2*nx*dx && pos.y() < 0.5*ny*dy);
  } else if (shape == "drop") {
    // a shallow pool of particles on the bottom
    double h = ny*dy/6.0;
    if (pos.y() < 2*h) return true;
    // and a sphere of particles above
    Vec3f center = Vec3f(nx*dx*0.5, 5*h,nz*dz*0.5);
    double length = (center-pos).Length();
    if (length < 0.8*h) return true;
    return false;
  } else {
    std::cout << "unknown shape: " << shape << std::endl;
    exit(0);
  }
}

// ==============================================================
void Fluid::GenerateParticles(const std::string &shape, const std::string &placement) {
  // create a set of points according to the "placement" token,
  // then check whether they are inside of the "shape"
  if (placement == "uniform") {
    int dens = (int)pow(density,0.334);
    assert (dens*dens*dens == density);
    // the uniform grid spacing
    double spacing = 1/double(dens);
    for (double x = 0.5*spacing*dx; x < nx*dx; x += spacing*dx) {
      for (double y = 0.5*spacing*dy; y < ny*dy; y += spacing*dy) {
        for (double z = 0.5*spacing*dz; z < nz*dz; z += spacing*dz) {
          Vec3f pos = Vec3f(x,y,z);
          if (inShape(pos,shape)) {
            Cell *cell = getCell(int(x/dx),int(y/dy),int(z/dz));
            FluidParticle *p = new FluidParticle();
            p->setPosition(pos);
            cell->addParticle(p);
          }
        }
      }
    }
  } else {
    assert (placement == "random");
    // note: we don't necessarily have the same number of particles in each cell
    for (int n = 0; n < nx*ny*nz*density; n++) {
      Vec3f pos = Vec3f(args->mtrand.rand()*nx*dx,
                        args->mtrand.rand()*ny*dy,
                        args->mtrand.rand()*nz*dz);
      if (inShape(pos,shape)) {      
        Cell *cell = getCell(int(pos.x()/dx),int(pos.y()/dy),int(pos.z()/dz));
        FluidParticle *p = new FluidParticle();
        p->setPosition(pos);
        cell->addParticle(p);
      }
    }
  }
}

// ==============================================================
// ==============================================================
// ANIMATION
// ==============================================================
// ==============================================================

void Fluid::Animate() {
  // MainMaster

  // the animation manager:  this is what gets done each timestep!

  // Updates the new_x_plus(), new values on back up. 
  // Old ones exist too
  ComputeNewVelocities();

  // This actually sets boundries
  SetBoundaryVelocities();
  
  // compressible / incompressible flow
  int count = 0;
  if (compressible == false) {
    double max_divergence = 0;
    for (int iters = 0; iters < 20; iters++) {
      count++;
      // What is diverance?
      if(Foster)
        max_divergence = AdjustForIncompressibility_Foster();
      else
        max_divergence = AdjustForIncompressibility();

      SetBoundaryVelocities();
      if (max_divergence < EPSILON) break;
    }
  }

  UpdatePressures();

  CopyVelocities();

  // advanced the particles through the fluid
  MoveParticles();
  ReassignParticles();

  // What does this do?
  SetEmptySurfaceFull();

  setupVBOs();
}

// ==============================================================

void Fluid::ComputeNewVelocities() {
  double dt = args->timestep;
  int i,j,k;

  // using the formulas from Foster & Metaxas
  // Updates to new_?_plus

  for (i = 0; i < nx-1; i++) {
    for (j = 0; j < ny; j++) {
      for (k = 0; k < nz; k++) {
        Cell *cell = getCell(i,j,k);
        double new_u_plus =
          get_u_plus(i,j,k) +            
          dt * ((1/dx) * (square(get_u_avg(i,j,k)) - square(get_u_avg(i+1,j,k))) +
                (1/dy) * (get_uv_plus(i,j-1,k) - get_uv_plus(i,j,k)) + 
                (1/dz) * (get_uw_plus(i,j,k-1) - get_uw_plus(i,j,k)) +
                args->gravity.x() +
                (1/dx) * (getPressure(i,j,k)-getPressure(i+1,j,k)) +
                (viscosity/square(dx)) * (get_u_plus(i+1,j  ,k  ) - 2*get_u_plus(i,j,k) + get_u_plus(i-1,j  ,k  )) +
                (viscosity/square(dy)) * (get_u_plus(i  ,j+1,k  ) - 2*get_u_plus(i,j,k) + get_u_plus(i  ,j-1,k  )) +
                (viscosity/square(dz)) * (get_u_plus(i  ,j  ,k+1) - 2*get_u_plus(i,j,k) + get_u_plus(i  ,j  ,k-1)) );
        cell->set_new_u_plus(new_u_plus);
      }
    }
  }

  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny-1; j++) {
      for (k = 0; k < nz; k++) {  
        Cell *cell = getCell(i,j,k);
        double new_v_plus =
          get_v_plus(i,j,k) +
          dt * ((1/dx) * (get_uv_plus(i-1,j,k) - get_uv_plus(i,j,k)) +
                (1/dy) * (square(get_v_avg(i,j,k)) - square(get_v_avg(i,j+1,k))) +
                (1/dz) * (get_vw_plus(i,j,k-1) - get_vw_plus(i,j,k)) +
                args->gravity.y() +
                (1/dy) * (getPressure(i,j,k)-getPressure(i,j+1,k)) +
                (viscosity/square(dx)) * (get_v_plus(i+1,j  ,k  ) - 2*get_v_plus(i,j,k) + get_v_plus(i-1,j  ,k  )) +
                (viscosity/square(dy)) * (get_v_plus(i  ,j+1,k  ) - 2*get_v_plus(i,j,k) + get_v_plus(i  ,j-1,k  )) +
                (viscosity/square(dz)) * (get_v_plus(i  ,j  ,k+1) - 2*get_v_plus(i,j,k) + get_v_plus(i  ,j  ,k-1)) );
        cell->set_new_v_plus(new_v_plus);
      }
    }
  }

  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++) {
      for (k = 0; k < nz-1; k++) {
        Cell *cell = getCell(i,j,k);
        double new_w_plus =
          get_w_plus(i,j,k) +
          dt * ((1/dx) * (get_uw_plus(i-1,j,k) - get_uw_plus(i,j,k)) +
                (1/dy) * (get_vw_plus(i,j-1,k) - get_vw_plus(i,j,k)) +
                (1/dz) * (square(get_w_avg(i,j,k)) - square(get_w_avg(i,j,k+1))) +
                args->gravity.z() +
                (1/dz) * (getPressure(i,j,k)-getPressure(i,j,k+1)) +
                (viscosity/square(dx)) * (get_w_plus(i+1,j  ,k  ) - 2*get_w_plus(i,j,k) + get_w_plus(i-1,j  ,k  )) +
                (viscosity/square(dy)) * (get_w_plus(i  ,j+1,k  ) - 2*get_w_plus(i,j,k) + get_w_plus(i  ,j-1,k  )) +
                (viscosity/square(dz)) * (get_w_plus(i  ,j  ,k+1) - 2*get_w_plus(i,j,k) + get_w_plus(i  ,j  ,k-1)) );
        cell->set_new_w_plus(new_w_plus);
      }
    }
  }
}


// ==============================================================

void Fluid::SetBoundaryVelocities() {

  // zero out flow perpendicular to the boundaries (no sources or sinks)
  for (int j = -1; j <= ny; j++) {
    for (int k = -1; k <= nz; k++) {
      getCell(-1  ,j,k)->set_u_plus(0);
      getCell(nx-1,j,k)->set_u_plus(0);
      getCell(nx  ,j,k)->set_u_plus(0);
    }
  }
  for (int i = -1; i <= nx; i++) {
    for (int k = -1; k <= nz; k++) {
      getCell(i,-1  ,k)->set_v_plus(0);
      getCell(i,ny-1,k)->set_v_plus(0);
      getCell(i,ny  ,k)->set_v_plus(0);
    }
  }
  for (int i = -1; i <= nx; i++) {
    for (int j = -1; j <= ny; j++) {
      getCell(i,j,-1  )->set_w_plus(0);
      getCell(i,j,nz-1)->set_w_plus(0);
      getCell(i,j,nz  )->set_w_plus(0);
    }
  }

  // free slip or no slip boundaries (friction with boundary)
  double xy_sign = (xy_free_slip) ? 1 : -1;
  double yz_sign = (yz_free_slip) ? 1 : -1;
  double zx_sign = (zx_free_slip) ? 1 : -1;
  for (int i = 0; i < nx; i++) {
    for (int j = -1; j <= ny; j++) {
      getCell(i,j,-1)->set_u_plus(xy_sign*getCell(i,j,0)->get_u_plus());
      getCell(i,j,nz)->set_u_plus(xy_sign*getCell(i,j,nz-1)->get_u_plus());
    }
    for (int k = -1; k <= nz; k++) {
      getCell(i,-1,k)->set_u_plus(zx_sign*getCell(i,0,k)->get_u_plus());
      getCell(i,ny,k)->set_u_plus(zx_sign*getCell(i,ny-1,k)->get_u_plus());
    }
  }
  for (int j = 0; j < ny; j++) {
    for (int i = -1; i <= nx; i++) {
      getCell(i,j,-1)->set_v_plus(xy_sign*getCell(i,j,0)->get_v_plus());
      getCell(i,j,nz)->set_v_plus(xy_sign*getCell(i,j,nz-1)->get_v_plus());
    }
    for (int k = -1; k <= nz; k++) {
      getCell(-1,j,k)->set_v_plus(yz_sign*getCell(0,j,k)->get_v_plus());
      getCell(nx,j,k)->set_v_plus(yz_sign*getCell(nx-1,j,k)->get_v_plus());
    }
  }
  for (int k = 0; k < nz; k++) {
    for (int i = -1; i <= nx; i++) {
      getCell(i,-1,k)->set_w_plus(zx_sign*getCell(i,0,k)->get_w_plus());
      getCell(i,ny,k)->set_w_plus(zx_sign*getCell(i,ny-1,k)->get_w_plus());
    }
    for (int j = -1; j <= ny; j++) {
      getCell(-1,j,k)->set_w_plus(yz_sign*getCell(0,j,k)->get_w_plus());
      getCell(nx,j,k)->set_w_plus(yz_sign*getCell(nx-1,j,k)->get_w_plus());
    }
  }
}

// ==============================================================

void Fluid::EmptyVelocities(int i, int j, int k) {
  Cell *c = getCell(i,j,k);
  if (c->getStatus() != CELL_EMPTY) return;
  Cell *ciplus = getCell(i+1,j,k);
  Cell *cjplus = getCell(i,j+1,k);
  Cell *ckplus = getCell(i,j,k+1);
  if (ciplus->getStatus() == CELL_EMPTY)
    c->set_new_u_plus(0);
  if (cjplus->getStatus() == CELL_EMPTY)
    c->set_new_v_plus(0);
  if (ckplus->getStatus() == CELL_EMPTY)
    c->set_new_w_plus(0);
}

void Fluid::CopyVelocities() {
  double dt = args->timestep;
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      for (int k = 0; k < nz; k++) {
  Cell *c = getCell(i,j,k);
  EmptyVelocities(i,j,k);
  c->copyVelocity();
  if (fabs(c->get_u_plus()) > 0.5*dx/dt ||
      fabs(c->get_v_plus()) > 0.5*dy/dt ||
      fabs(c->get_w_plus()) > 0.5*dz/dt) {
    // velocity has exceeded reasonable threshhold
    std::cout << "velocity has exceeded reasonable threshhold, stopping animation" << std::endl;
    args->animate=false;
  }
      }
    }
  }
}

// ==============================================================
double Fluid::AdjustForIncompressibility_Foster(){
  // Unlike code in class we use approximation methods described in the
  // foster paper to calcualte the divergance taking into account beta,
  // a relaxation coffienent. We Also this this bindly for all cells,
  // while zeroing out the empty cells after.
  // note it takes 2-5 iterations to converage on some episone descripbed
  // the paper.


 double max_divergence = -1;

  for (int i = -1; i <= nx; i++) {
    for (int j = -1; j <= ny; j++) {
      for (int k = -1; k <= nz; k++) {
        // For each cell in the system

        if (i >= 0 && i < nx && j >= 0 && j < ny && k >= 0 && k < nz) {

          // Get divergence
          double divergence = 
          - ( (1/dx) * (get_new_u_plus(i,j,k) - get_new_u_plus(i-1,j,k)) +
          (1/dy) * (get_new_v_plus(i,j,k) - get_new_v_plus(i,j-1,k)) +
          (1/dz) * (get_new_w_plus(i,j,k) - get_new_w_plus(i,j,k-1)) );

          if(max_divergence < fabs(divergence)) max_divergence = fabs(divergence);

          // Take the change in time into concideration
          double dt = args->timestep;

          // Our relaxation coffient taken into account
          double beta = BETA_0/((2*dt) * (1/square(dx) + 1/square(dy) + 1/square(dz)));
          double dp = beta*(divergence);

          // Update the velocities (based on the equations 5-7 provdied in paper)
          adjust_new_u_plus(i,j,k, (dt/dx)*dp);
          adjust_new_u_plus(i-1,j,k, -1*(dt/dx)*dp);

          adjust_new_v_plus(i,j,k, (dt/dy)*dp);
          adjust_new_v_plus(i,j-1,k, -1*(dt/dy)*dp);

          adjust_new_w_plus(i,j,k, (dt/dz)*dp);
          adjust_new_w_plus(i,j,k-1, -1*(dt/dz)*dp);



        }else{
          // We have a boundry cell




        }
      }

    }//for each cell
  }

  return max_divergence;
}



double Fluid::AdjustForIncompressibility() {
  // JUMP

 double max_divergence = -1;

  for (int i = -1; i <= nx; i++) {
    for (int j = -1; j <= ny; j++) {
      for (int k = -1; k <= nz; k++) {
        // For each cell in the system
        Cell *c = getCell(i,j,k);

        if (i >= 0 && i < nx && j >= 0 && j < ny && k >= 0 && k < nz) {

          // If you're not full skip
          if(c->getStatus() != CELL_FULL)
            continue;

          // Get divergence
          double divergence = 
          - ( (1/dx) * (get_new_u_plus(i,j,k) - get_new_u_plus(i-1,j,k)) +
          (1/dy) * (get_new_v_plus(i,j,k) - get_new_v_plus(i,j-1,k)) +
          (1/dz) * (get_new_w_plus(i,j,k) - get_new_w_plus(i,j,k-1)) );

          if(max_divergence < fabs(divergence)) max_divergence = fabs(divergence);

          //Spread the unhappyness to "full" cells
          double chunk = divergence / (double) getLegalAdjCells(i,j,k);

          //Stats /////////////////////
          // printf("-=-=-=-=-=--=-=\n");
          // printf("-=-=-[%d,%d,%d]-=-=-=-\n",i,j,k);
          // printf("divergence:  %f\n",divergence);
          // printf("divergance chunk: %f\n", chunk);

          // East Face
          if(legal_full_cell(i+1,j,k))
            adjust_new_u_plus(i,j,k,chunk);

          // South Face
          if(legal_full_cell(i,j+1,k))
            adjust_new_v_plus(i,j,k,chunk);

          // Face-me
          if(legal_full_cell(i,j,k+1))
            adjust_new_w_plus(i,j,k,chunk);

          // West Face
          if(legal_full_cell(i-1,j,k))
            adjust_new_u_plus(i-1,j,k, -1*chunk);

          // North Face
          if(legal_full_cell(i,j-1,k))
            adjust_new_v_plus(i,j-1,k, -1*chunk);

          // Face-Away
          if(legal_full_cell(i,j,k-1))
            adjust_new_w_plus(i,j,k-1,-1*chunk);
        }
      }
    }
  }
  return max_divergence;
}

// ==============================================================

void Fluid::UpdatePressures() {
  for (int i = -1; i <= nx; i++) {
    for (int j = -1; j <= ny; j++) {
      for (int k = -1; k <= nz; k++) {
        // For each cell in the system
        Cell *c = getCell(i,j,k);
        if (i >= 0 && i < nx && j >= 0 && j < ny && k >= 0 && k < nz) {
          // compute divergence and increment/decrement pressure
          double pressure = c->getPressure();
          double divergence = 
            - ( (1/dx) * (get_new_u_plus(i,j,k) - get_new_u_plus(i-1,j,k)) +
          (1/dy) * (get_new_v_plus(i,j,k) - get_new_v_plus(i,j-1,k)) +
          (1/dz) * (get_new_w_plus(i,j,k) - get_new_w_plus(i,j,k-1)) );
          double dt = args->timestep;
          double beta = BETA_0/((2*dt) * (1/square(dx) + 1/square(dy) + 1/square(dz)));
          double dp = beta*divergence;
          c->setPressure(pressure + dp);
        } else {
          // zero out boundary cells (just in case)
          c->setPressure(0);
        }

        // zero out empty cells (From Foster 2001 paper)
        if (c->getStatus() == CELL_EMPTY) {
          c->setPressure(0);
        }

        // ========================================

      }
    }
  }
}

// ==============================================================

void Fluid::MoveParticles() {

  // Get timestep 
  double dt = args->timestep;

  // For each particle
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      for (int k = 0; k < nz; k++) {

        Cell *cell = getCell(i,j,k);

        // Get particles inside cell
        std::vector<FluidParticle*> &particles = cell->getParticles();

        for (unsigned int iter = 0; iter < particles.size(); iter++) {

          // For each particle in cell
          FluidParticle *p = particles[iter];
          Vec3f pos = p->getPosition();
          Vec3f vel = getInterpolatedVelocity(pos);

          // Update
          Vec3f pos2 = pos + vel*dt;


          // euler integration
          p->setPosition(pos2);

          // I will halve the timestemp if i realize we are moving to fast, so the next iteration won't be as effected
          if(fabs(pos.Distance3f(pos2)) >= fabs(dx)){
            args->timestep = args->timestep/2;
            std::cout << "Happening to fast, Time Halved -->" << args->timestep << std::endl;
          }

        }
      }
    }
  }
}

// ==============================================================

void Fluid::ReassignParticles() {
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      for (int k = 0; k < nz; k++) {
        Cell *cell = getCell(i,j,k);
  std::vector<FluidParticle*> &particles = cell->getParticles();
        for (unsigned int iter = 0; iter < particles.size(); iter++) {
          FluidParticle *p = particles[iter];
          Vec3f pos = p->getPosition();
          int i2 = (int)my_min(double(nx-1),my_max(0.0,floor(pos.x()/dx)));
          int j2 = (int)my_min(double(ny-1),my_max(0.0,floor(pos.y()/dy)));
          int k2 = (int)my_min(double(nz-1),my_max(0.0,floor(pos.z()/dz)));
          // if the particle has crossed one of the cell faces 
          // assign it to the new cell
          if (i != i2 || j != j2 || k != k2) {
            cell->removeParticle(p);
            getCell(i2,j2,k2)->addParticle(p);
          } 
        }
      }
    }
  }
}

// ==============================================================

void Fluid::SetEmptySurfaceFull() {
  int i,j,k;
  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++) {
      for (k = 0; k < nz; k++) {
        Cell *cell = getCell(i,j,k);
        if (cell->numParticles() == 0)
          cell->setStatus(CELL_EMPTY);
        else 
          cell->setStatus(CELL_FULL);
      }
    }
  }

  // pick out the boundary cells
  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++) {
      for (k = 0; k < nz; k++) {
        Cell *cell = getCell(i,j,k);
        if (cell->getStatus() == CELL_FULL &&
            (getCell(i-1,j,k)->getStatus() == CELL_EMPTY ||
             getCell(i+1,j,k)->getStatus() == CELL_EMPTY ||
             getCell(i,j-1,k)->getStatus() == CELL_EMPTY ||
             getCell(i,j+1,k)->getStatus() == CELL_EMPTY ||
             getCell(i,j,k-1)->getStatus() == CELL_EMPTY ||
             getCell(i,j,k+1)->getStatus() == CELL_EMPTY)) {
          cell->setStatus(CELL_SURFACE);
        }
      }
    }
  }
}

// ==============================================================

Vec3f Fluid::getInterpolatedVelocity(const Vec3f &pos) const {
  // Input: given a single postion in our space
  // Output: the interpola2: b = (Vec3f &) @0x7fffffffc010: {data = {0, 45, 0}}

  // TODO

  int i = int(floor(pos.x()/dx)); if (i < 0) i = 0; if (i >= nx) i = nx-1;
  int j = int(floor(pos.y()/dy)); if (j < 0) j = 0; if (j >= ny) j = ny-1;
  int k = int(floor(pos.z()/dz)); if (k < 0) k = 0; if (k >= nz) k = nz-1;

  //printf("[%.2f,%.2f,%.2f] Particle Interpolation Vel____________\n",pos.x(),pos.y(),pos.z());

  // What points are we going to use?
  const Vec3f corner(i*dx, j*dy, k*dz);
  //printf("[%.2f,%.2f,%.2f] Assigned Corner\n",i*dx,j*dy,k*dz);
  // Horizontal Velcoities (u) //////////////////////////////////////////

  Vec3f pos_u1,pos_u2,pos_u3,pos_u4;
  double vel_u1,vel_u2,vel_u3,vel_u4;

  // Get first 2
  pos_u1 = {corner.x(),      corner.y() + (0.5 * dy), corner.z()};
  pos_u2 = {corner.x() + dx, corner.y() + (0.5 * dy), corner.z()};

  vel_u1 = get_u_plus(i-1,j,k);
  vel_u2 = get_u_plus(i,j,k);

  // find the remaining last two data sets

  if( pos.y() > corner.y() + (0.5* dy) ){
    // If i lie below a threshold
    pos_u3 = {corner.x(),      corner.y() + dy + (0.5 * dy), corner.z() };
    pos_u4 = {corner.x() + dx, corner.y() + dy + (0.5 * dy), corner.z() };
    
    vel_u3 = get_u_plus(i-1, j+1, k);
    vel_u4 = get_u_plus(i, j + 1, k);

  }else{  

    // If i lie above
    pos_u3 = {corner.x(),      corner.y() - (0.5 * dy), corner.z() };
    pos_u4 = {corner.x() + dx, corner.y() - (0.5 * dy), corner.z() };
    
    vel_u3 = get_u_plus(i-1,j-1,k);
    vel_u4 = get_u_plus(i,  j-1,k);

  }

  // printf("Horz Influ [%.2f,%.2f] is %.2f\n", pos_u1.x(),pos_u1.y(), getAreaSquares(pos,pos_u1));
  // printf("Horz Influ [%.2f,%.2f] is %.2f\n", pos_u2.x(),pos_u2.y(), getAreaSquares(pos,pos_u2));
  // printf("Horz Influ [%.2f,%.2f] is %.2f\n", pos_u3.x(),pos_u3.y(), getAreaSquares(pos,pos_u3));
  // printf("Horz Influ [%.2f,%.2f] is %.2f\n", pos_u4.x(),pos_u4.y(), getAreaSquares(pos,pos_u4));
  //printf("Total Horz Area %f\n", getAreaSquares(pos,pos_u1) + getAreaSquares(pos,pos_u2) + getAreaSquares(pos,pos_u3) + getAreaSquares(pos,pos_u4));
  //printf("Goal  Horz Area %f\n", dx*dy);

  double interpolated_u = getAreaSquares(pos, pos_u1) * vel_u1;
  interpolated_u       += getAreaSquares(pos, pos_u2) * vel_u2;
  interpolated_u       += getAreaSquares(pos, pos_u3) * vel_u3;
  interpolated_u       += getAreaSquares(pos, pos_u4) * vel_u4;
  interpolated_u       *= (1/(dx*dy));


  // Verteical Velcoities (v) //////////////////////////////////////////
  Vec3f pos_v1,pos_v2,pos_v3,pos_v4;
  double vel_v1,vel_v2,vel_v3,vel_v4;


  // Get first 2
  pos_v1 = {corner.x() + (0.5 * dx), corner.y(),      corner.z()};
  pos_v2 = {corner.x() + (0.5 * dx), corner.y() + dy, corner.z()};

  vel_v1 = get_v_plus(i, j-1, k);
  vel_v2 = get_v_plus(i, j, k );

  // find the remaining last two data sets
  if( pos.x() > corner.x() + (0.5* dx) ){

    // If i lie above
    pos_v3 = {corner.x() + dx + (0.5 * dx), corner.y(),      corner.z() };
    pos_v4 = {corner.x() + dx + (0.5 * dx), corner.y() + dy, corner.z() };
    
    vel_v3 = get_v_plus(i+1, j-1,   k);
    vel_v4 = get_v_plus(i+1, j, k);

  }else{
    // If i lie below a threshold
    pos_v3 = {corner.x() - dx, corner.y(),      corner.z() };
    pos_v4 = {corner.x() - dx, corner.y() + dy, corner.z() };
    
    vel_v3 = get_v_plus(i-1, j-1,   k);
    vel_v4 = get_v_plus(i-1, j, k);

  }

  //printf("-------------------------------\n");

  //printf("Total Vert Area %f\n", getAreaSquares(pos,pos_v1) + getAreaSquares(pos,pos_v2) + getAreaSquares(pos,pos_v3) + getAreaSquares(pos,pos_v4));
  //printf("Goal  Vert Area %f\n", dx*dy);
  // printf("Vert Influ [%.2f%.2f] is %.2f\n", pos_v1.x(), pos_v1.y(), getAreaSquares(pos,pos_v1));
  // printf("Vert Influ [%.2f,%.2f] is %.2f\n", pos_v2.x(), pos_v2.y(), getAreaSquares(pos,pos_v2));
  // printf("Vert Influ [%.2f,%.2f] is %.2f\n", pos_v3.x(), pos_v3.y(), getAreaSquares(pos,pos_v3));
  // printf("Vert Influ [%.2f,%.2f] is %.2f\n", pos_v4.x(), pos_v4.y(), getAreaSquares(pos,pos_v4));

  double interpolated_v = getAreaSquares(pos, pos_v1) * vel_v1;
  interpolated_v       += getAreaSquares(pos, pos_v2) * vel_v2;
  interpolated_v       += getAreaSquares(pos, pos_v3) * vel_v3;
  interpolated_v       += getAreaSquares(pos, pos_v4) * vel_v4;
  interpolated_v       *= (1/(dx*dy));

  // 3D Velcoities (?) //////////////////////////////////////////

  return Vec3f(interpolated_u,interpolated_v,get_w_avg(i,j,k));
  // *********************************************************************  
}


double Fluid::getAreaSquares(const Vec3f& a, Vec3f&b) const {
  // Assumption Vec3f is 

  Vec3f low_bound_a(a.x() - dx/2, a.y() - dy/2, a.z() - dz/2);
  Vec3f hig_bound_a(a.x() + dx/2, a.y() + dy/2, a.z() + dz/2);

  Vec3f low_bound_b(b.x() - dx/2, b.y() - dy/2, b.z() - dz/2);
  Vec3f hig_bound_b(b.x() + dx/2, b.y() + dy/2, b.z() + dz/2);

  double width = 0;
  double length = 0;

  // X's intersect
  if(low_bound_a.x() < low_bound_b.x() && low_bound_b.x() <= hig_bound_a.x()){
    // a is first
    width = hig_bound_a.x() - low_bound_b.x();

  }else if(low_bound_a.x() < hig_bound_b.x() && hig_bound_b.x() <= hig_bound_a.x()){
    // b is first
    width = hig_bound_b.x() - low_bound_a.x();

  }else if( low_bound_a.x() == low_bound_b.x() && hig_bound_a.x() == hig_bound_b.x() ){
    // a and b are the same, overlap
    width = dx;
  }else{
    width = 0;
  }

  // Y's intersect
  if(low_bound_a.y() < low_bound_b.y() && low_bound_b.y() <= hig_bound_a.y()){
    // a is first
    length = hig_bound_a.y() - low_bound_b.y();

  }else if(low_bound_a.y() < hig_bound_b.y() && hig_bound_b.y() <= hig_bound_a.y()){
    // b is first
    length = hig_bound_b.y() - low_bound_a.y();

  }else if( low_bound_a.y() == low_bound_b.y() && hig_bound_a.y() == hig_bound_b.y() ){
    // a and b are the same, overlap
    length = dy;
  }else{
    length=0;
  }

  assert(width*length <= dx*dy);
  assert(width*length >= 0);
  return width*length;

}


