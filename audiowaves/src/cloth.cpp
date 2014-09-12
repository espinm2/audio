#include "glCanvas.h"
#include <fstream>
#include <cmath>
#include <math.h>       /* fabs */>
#include "cloth.h"
#include "argparser.h"
#include "vectors.h"
#include "utils.h"
using std::vector;

// ================================================================================
// TODO: Implement improved timestep fix
//       * Define a the limit as a weird jump in points after the first iteration.
//       * Half the time step
// TODO: Create new cloth case
//       * Large cloth with the corner fixed, let hung but one
// ================================================================================

Cloth::Cloth(ArgParser *_args) {
  args =_args;

  // open the file
  std::ifstream istr(args->cloth_file.c_str());
  assert (istr != NULL);
  std::string token;

  // read in the simulation parameters
  istr >> token >> k_structural; assert (token == "k_structural");  // (units == N/m)  (N = kg*m/s^2)
  istr >> token >> k_shear; assert (token == "k_shear");
  istr >> token >> k_bend; assert (token == "k_bend");
  istr >> token >> damping; assert (token == "damping");
  // NOTE: correction factor == .1, means springs shouldn't stretch more than 10%
  //       correction factor == 100, means don't do any correction
  istr >> token >> provot_structural_correction; assert (token == "provot_structural_correction");
  istr >> token >> provot_shear_correction; assert (token == "provot_shear_correction");

  // the cloth dimensions
  istr >> token >> nx >> ny; // (units == meters)
  assert (token == "m");
  assert (nx >= 2 && ny >= 2);

  // the corners of the cloth
  Vec3f a,b,c,d;
  istr >> token >> a; assert (token == "p");
  istr >> token >> b; assert (token == "p");
  istr >> token >> c; assert (token == "p");
  istr >> token >> d; assert (token == "p");

  // fabric weight  (units == kg/m^2)
  // denim ~300 g/m^2
  // silk ~70 g/m^2
  double fabric_weight;
  istr >> token >> fabric_weight; assert (token == "fabric_weight");
  double area = AreaOfTriangle(a,b,c) + AreaOfTriangle(a,c,d);

  // create the particles
  particles = new ClothParticle[nx*ny];
  double mass = area*fabric_weight / double(nx*ny);
  for (int i = 0; i < nx; i++) {
    double x = i/double(nx-1);
    Vec3f ab = (1-x)*a + x*b;
    Vec3f dc = (1-x)*d + x*c;
    for (int j = 0; j < ny; j++) {
      double y = j/double(ny-1);
      ClothParticle &p = getParticle(i,j);
      Vec3f abdc = (1-y)*ab + y*dc;
      p.setOriginalPosition(abdc);
      p.setPosition(abdc);
      p.setVelocity(Vec3f(0,0,0));
      p.setMass(mass);
      p.setFixed(false);
    }
  }

  // the fixed particles
  while (istr >> token) {
    assert (token == "f");
    int i,j;
    double x,y,z;
    istr >> i >> j >> x >> y >> z;
    ClothParticle &p = getParticle(i,j);
    p.setPosition(Vec3f(x,y,z));
    p.setFixed(true);
  }

  computeBoundingBox();
  initializeVBOs();
  setupVBOs();
}

// ================================================================================

void Cloth::computeBoundingBox() {
  box = BoundingBox(getParticle(0,0).getPosition());
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      box.Extend(getParticle(i,j).getPosition());
      box.Extend(getParticle(i,j).getOriginalPosition());
    }
  }
}

// ================================================================================

const vector<ClothParticle*> Cloth::getAdjParticles(int i, int j) {
  // Input: Given an i,j index
  // Output: This will return the adj (structural) particles in a vector

  //Check if valid index
  assert(0 <= i && i < nx);
  assert(0 <= j && j < nx);

  vector<ClothParticle*> adjPartVec;
  ClothParticle * ptr = NULL; 

  // Get left
  if(j - 1 >= 0){
    ptr = &getParticle(i,j-1);
    adjPartVec.push_back(ptr);
  }
  // Get right
  if(j + 1 < ny){
    ptr = &getParticle(i,j+1);
    adjPartVec.push_back(ptr);
  }
  // Get top
  if(i - 1 >= 0){
    ptr = &getParticle(i-1,j);
    adjPartVec.push_back(ptr);
  }
  // Get bottom
  if(i + 1 < nx){
    ptr = &getParticle(i+1,j);
    adjPartVec.push_back(ptr);
  }

  return  adjPartVec;
}

const vector<ClothParticle*> Cloth::getShearParticles(int i, int j) {
  // Input: Given an i,j index
  // Output: This will return the adj (Shear) particles in a vector

  //Check if valid index
  assert(0 <= i && i < nx);
  assert(0 <= j && j < nx);

  vector<ClothParticle*> shearPartVec;
  ClothParticle * ptr = NULL; 

  // Get upperLeft
  if(i - 1 >= 0 && j - 1 >= 0){
    ptr = &getParticle( i-1 , j-1);
    shearPartVec.push_back(ptr);
  }

  // Get lowerRight
  if(i + 1 < nx && j + 1 < ny){
    ptr = &getParticle( i + 1 , j + 1);
    shearPartVec.push_back(ptr);
  }

  // Get upperRight
  if(i - 1 >= 0 && j + 1 < ny){
    ptr = &getParticle(i-1, j+1);
    shearPartVec.push_back(ptr);
  }

  // Get lowerLeft
  if(i + 1 < nx && j - 1 >= 0){
    ptr = &getParticle(i+1,j-1);
    shearPartVec.push_back(ptr);
  }

  return shearPartVec;
}

// ================================================================================

const vector<ClothParticle*> Cloth::getFlexParticles(int i, int j) {

  //Check if valid index
  assert(0 <= i && i < nx);
  assert(0 <= j && j < nx);

  vector<ClothParticle*> flexPartVec;
  ClothParticle * ptr = NULL; 

  // Get top
  if(i - 2 >= 0){
    ptr = &getParticle(i-2,j);
    flexPartVec.push_back(ptr);
  }

  // Get bot
  if(i + 2 < nx){
    ptr = &getParticle(i+2,j);
    flexPartVec.push_back(ptr);
  }


  // Get left
  if(j - 2 >= 0){
    ptr = &getParticle( i  , j - 2);
    flexPartVec.push_back(ptr);
  }

  // Get right
  if(j + 2 < ny ){
    ptr = &getParticle( i  , j + 2);
    flexPartVec.push_back(ptr);
  }

  return flexPartVec;
}

// ================================================================================

const Vec3f Cloth::getSpringForce(ClothParticle* a, ClothParticle* b) {

  //TODO Calc spring const between these two
  Vec3f aOrginalPos = a->getOriginalPosition();
  Vec3f bOrginalPos = b->getOriginalPosition();
  double restLength = aOrginalPos.Distance3f(bOrginalPos);
  Vec3f p_i = a->getPosition();
  Vec3f p_j = b->getPosition();

  // JUMP
  double displace = k_structural * (p_i.Distance3f(p_j) - restLength);
  Vec3f ratio = (p_j - p_i) * (1/p_i.Distance3f(p_j));
  return displace * ratio;
}

// ================================================================================

void Cloth::Animate() {

  // Get the fixed distance for  structure particles
  double restStructLength = getParticle(0,0).getOriginalPosition().Distance3f(getParticle(0,1).getOriginalPosition());
  double restShearLength = getParticle(0,0).getOriginalPosition().Distance3f(getParticle(1,1).getOriginalPosition());

  // For each particle in the system update the states
  for( int i = 0; i < nx; i++){

    for( int j = 0; j < ny; j++){

      // Trying to do eulers method
      ClothParticle* curP = &getParticle(i,j);

      // Unless fixed
      if((curP->isFixed()))
        continue;

      // Calculating Forces///////////////////////////////
      // * Gravity
      // * Spring - Structural,Shear,Bend
      // * Dampening 

      Vec3f f_damp = -1* damping * curP->getVelocity();
      Vec3f f_gravity = curP->getMass() * args->gravity;
      // f_spring = addition of all forces
      Vec3f f_spring;

      // Structural
      vector<ClothParticle*> adjPartVec = getAdjParticles(i,j);
      for(unsigned int v = 0; v < adjPartVec.size(); v++)
        f_spring = f_spring + getSpringForce(curP,adjPartVec[v]);
      // Shear 
      vector<ClothParticle*> shearVec = getShearParticles(i,j);
      for(unsigned int v = 0; v < shearVec.size(); v++)
        f_spring = f_spring + getSpringForce(curP,shearVec[v]);
      // Flex
      vector<ClothParticle*> flexVec = getFlexParticles(i,j);
      for(unsigned int v = 0; v < flexVec.size(); v++)
        f_spring = f_spring + getSpringForce(curP,flexVec[v]);

      // Accleration
      curP->setAcceleration((1/curP->getMass()) * (f_spring + f_gravity + f_damp ));

      // Update Velocity
      curP->setVelocity(curP->getVelocity() + (args->timestep * curP->getAcceleration()));

      // Get last Position to check against
      Vec3f oldPos = curP->getPosition();
      // Update Position
      curP->setPosition(curP->getPosition() + args->timestep * curP->getVelocity());


      // I will halve the timestemp if i realize we are moving to fast, so the next iteration won't be as effected
      if(fabs(oldPos.Distance3f(curP->getPosition())) >= fabs(restStructLength/4.0)){
        //std::cout << "Happening to fast, Time Halved\n";
        args->timestep = args->timestep/2;
      }

    }
  }


  // For each particle in the system check if edges around it is overstreched
  for(unsigned int i = 0; i < nx; i++){

    for(unsigned int j = 0; j < ny; j++){

      // Collecting surrounding nodes 
      ClothParticle* curP = &getParticle(i,j);
      vector<ClothParticle*> adjPartVec = getAdjParticles(i,j);

      for(unsigned int v = 0; v < adjPartVec.size(); v++){
        // For edge check if oversteched
        ClothParticle * otherP = adjPartVec[v];

        double distance = std::abs(curP->getPosition().Distance3f(otherP->getPosition()));
        double max_struct_length = (1 + provot_structural_correction) *restStructLength;
        double diff_distance = std::abs(max_struct_length - distance);

        if(max_struct_length - distance < 0 && provot_structural_correction != 100){
          // std::cout << "-=-=-=-=-= Struct Correction -=-=-=-=-=-=-==-=\n";
          // std::cout << "Overstreched at Index: (" << i << "," << j << ")\n";
          // std::cout << "Spring length: " << distance << std::endl;
          // std::cout << "Spring MaxLen: " << max_struct_length << std::endl;
          // std::cout << "Correction  : " << diff_distance << std::endl;

          Vec3f line_cur_other = otherP->getPosition() - curP->getPosition();
          Vec3f line_other_cur = curP->getPosition() - otherP->getPosition();
          Vec3f new_cur = curP->getPosition();
          Vec3f new_other = otherP->getPosition();


          // Which are fixed? none, one or both?
          if(curP->isLoose() && otherP->isLoose()){
            // std::cout << "Both Loose\n";
            // Move both
            new_cur = curP->getPosition() + ((diff_distance/2.0) * ((1/line_cur_other.Length()) * line_cur_other));
            new_other = otherP->getPosition() +((diff_distance/2.0) *((1/line_other_cur.Length()) * line_other_cur));

            // std::cout << "Target:        " << max_struct_length << std::endl;
            // std::cout << "Result:        " << std::abs(new_cur.Distance3f(new_other))<<std::endl;

          }else if( curP->isLoose() && otherP->isFixed()){
            // std::cout << "One Loose\n";


            new_cur = curP->getPosition() + ((diff_distance) *((1/line_cur_other.Length()) * line_cur_other));
            // std::cout << "Target:        " << max_struct_length << std::endl;
            // std::cout << "Result:        " << std::abs(new_cur.Distance3f(new_other))<<std::endl;

          }else if( curP->isFixed() && otherP->isLoose()){
            // std::cout << "One Loose\n";

            new_other = otherP->getPosition() +((diff_distance) *((1/line_other_cur.Length()) * line_other_cur));
            // std::cout << "Target:        " << max_struct_length << std::endl;
            // std::cout << "Result:        " << std::abs(new_cur.Distance3f(new_other)) << std::endl;

          }else{
            // Then both are fixed, in which case do nothing.
            assert(curP->isFixed() && otherP->isFixed());
          }

          curP->setPosition(new_cur);
          otherP->setPosition(new_other);
        }
      }



      //SHEAR
      vector<ClothParticle*> shearPartVec = getShearParticles(i,j);
      for(unsigned int v = 0; v < shearPartVec.size(); v++){
        // For edge check if oversteched
        ClothParticle * otherP = shearPartVec[v];

        double distance = std::abs(curP->getPosition().Distance3f(otherP->getPosition()));
        double max_shear_length = (1+ provot_shear_correction) *restShearLength;
        double diff_distance = std::abs(max_shear_length - distance);


        if(max_shear_length - distance < 0 && provot_shear_correction != 100){

          // std::cout << "-=-=-=-=-= Shear Correction -=-=-=-=-=-=-==-=\n";
          // std::cout << "Overstreched at Index: (" << i << "," << j << ")\n";
          // std::cout << "Spring length: " << distance << std::endl;
          // std::cout << "Spring MaxLen: " << max_shear_length << std::endl;
          // std::cout << "Correction   : " << diff_distance << std::endl;

          Vec3f line_cur_other = otherP->getPosition() - curP->getPosition();
          Vec3f line_other_cur = curP->getPosition() - otherP->getPosition();
          Vec3f new_cur = curP->getPosition();
          Vec3f new_other = otherP->getPosition();


          // Which are fixed? none, one or both?
          if(curP->isLoose() && otherP->isLoose()){
            // Move both
            new_cur = curP->getPosition() + ((diff_distance/2.0) * ((1/line_cur_other.Length()) * line_cur_other));
            new_other = otherP->getPosition() +((diff_distance/2.0) *((1/line_other_cur.Length()) * line_other_cur));
            // std::cout << "Target:        " << max_shear_length << std::endl;
            // std::cout << "Result:        " << std::abs(new_cur.Distance3f(new_other))<<std::endl;


          }else if( curP->isLoose() && otherP->isFixed()){
            new_cur = curP->getPosition() + ((diff_distance) *((1/line_cur_other.Length()) * line_cur_other));
            // std::cout << "Target:        " << max_shear_length << std::endl;
            // std::cout << "Result:        " << std::abs(new_cur.Distance3f(new_other))<<std::endl;

          }else if( curP->isFixed() && otherP->isLoose()){
            new_other = otherP->getPosition() +((diff_distance) *((1/line_other_cur.Length()) * line_other_cur));
            // std::cout << "Target:        " << max_shear_length << std::endl;
            // std::cout << "Result:        " << std::abs(new_cur.Distance3f(new_other))<<std::endl;

          }else{
            // Then both are fixed, in which case do nothing.
            assert(curP->isFixed() && otherP->isFixed());
          }

          curP->setPosition(new_cur);
          otherP->setPosition(new_other);
        }
      }
    }
  }
  // redo VBOs for rendering
  setupVBOs();
}

// ================================================================================
