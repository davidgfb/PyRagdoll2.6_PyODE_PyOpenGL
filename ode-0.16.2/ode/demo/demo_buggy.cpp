/*
buggy with suspension.
this also shows you how to use geom groups.
*/

#include <ode/ode.h>
#include <drawstuff/drawstuff.h>
#include "texturepath.h"

#ifdef _MSC_VER
	#pragma warning(disable:4244 4305)  // for VC++, no precision loss complaints
#endif

// select correct drawing functions
#ifdef dDOUBLE
	#define dsDrawBox dsDrawBoxD
	#define dsDrawSphere dsDrawSphereD
	#define dsDrawCylinder dsDrawCylinderD
	#define dsDrawCapsule dsDrawCapsuleD
#endif

// some constants
#define LENGTH 0.7	// chassis length
#define WIDTH 0.5	// chassis width
#define HEIGHT 0.2	// chassis height
#define RADIUS 0.18	// wheel radius
#define STARTZ 0.5	// starting height of chassis
#define CMASS 1		// chassis mass
#define WMASS 0.2	// wheel mass

static const dVector3 yunit = { 0, 1, 0 }, zunit = { 0, 0, 1 };


// dynamics and collision objects (chassis, 3 wheels, environment)
static dWorldID world;
static dSpaceID space;
static dBodyID body[4];
static dJointID joint[3];	// joint[0] is the front wheel
static dJointGroupID contactgroup;
static dGeomID ground;
static dSpaceID car_space;
static dGeomID box[1];
static dGeomID sphere[3];
static dGeomID ground_box;


// things that the user controls
static dReal speed = 0, steer = 0;	// user commands


// this is called by dSpaceCollide when two objects in space are
// potentially colliding.
static void nearCallback (void *, dGeomID o1, dGeomID o2)
{
  int i, n;

  // only collide things with the ground
  int g1 = (o1 == ground || o1 == ground_box); // booleanos
  int g2 = (o2 == ground || o2 == ground_box);
  
  if (!(g1 ^ g2)) {
  	return;
  }

  const int N = 10;
  dContact contact[N];
  n = dCollide (o1,o2,N,&contact[0].geom,sizeof(dContact));
  
  if (n > 0) {
    for (i=0; i<n; i++) {
      contact[i].surface.mode = dContactSlip1 | dContactSlip2 |
	dContactSoftERP | dContactSoftCFM | dContactApprox1;
      
      contact[i].surface.mu = dInfinity;
      
      contact[i].surface.slip1    = 0.1;
      contact[i].surface.slip2    = 0.1;
      contact[i].surface.soft_erp = 0.5;
      contact[i].surface.soft_cfm = 0.3;
      
      dJointID c = dJointCreateContact (world,contactgroup,&contact[i]);
      dJointAttach (c,
		    dGeomGetBody(contact[i].geom.g1),
		    dGeomGetBody(contact[i].geom.g2));
    }
  }
}


// start simulation - set viewpoint

static void start()
{
  dAllocateODEDataForThread(dAllocateMaskAll);

  static float xyz[3] = {0.8317f, -0.9817f, 0.8f};
  static float hpr[3] = { 121.0f,   -27.5f, 0.0f};
  dsSetViewpoint (xyz,hpr);
  printf ("Press:\t'a' to increase speed.\n"
	  "\t'z' to decrease speed.\n"
	  "\t',' to steer left.\n"
	  "\t'.' to steer right.\n"
	  "\t' ' to reset speed and steering.\n"
	  "\t'1' to save the current state to 'state.dif'.\n");
}


// called when a key pressed //TODO arregla limites giro
static void command (int cmd)
{
	if (cmd == 'w' || cmd == 'W') { // lower // flechas 
		speed += 0.3;
	}
	
	if (cmd == 'a' || cmd == 'A') {
		if (steer > -1.0) {
	    		steer -= 0.5;
    		}
    	
    		//printf("%.6f\n", steer);
	}
	
	if (cmd == 's' || cmd == 'S') {
		speed -= 0.3;
	}
	
	if (cmd == 'd' || cmd == 'D') {
		if (steer < 1.0) {
    			steer += 0.5;
    		}
    	
    		//printf("%.6f\n", steer);
	}
	
  	switch (cmd) { //lower      
  		case ' ': // espacio = freno de mano
    			speed = 0; //speed, steer = 0, 0; no funciona tampoco (0, 0)
    			steer = 0;
    		break;
  	}
}

// simulation loop
static void simLoop (int pause)
{
  int i;
  
  if (!pause) {
    // motor
    dJointSetHinge2Param (joint[0],  dParamVel2, -speed);
    dJointSetHinge2Param (joint[0], dParamFMax2,    0.1);

    // steering
    dReal v = steer - dJointGetHinge2Angle1 (joint[0]);
    
    if (v > 0.1) {
    	v = 0.1;
    }
    
    if (v < -0.1) { 
    	v = -0.1;
    }
    
    v *= 10.0;
    dJointSetHinge2Param (joint[0], dParamVel,           v);
    dJointSetHinge2Param (joint[0], dParamFMax,        0.2);
    dJointSetHinge2Param (joint[0], dParamLoStop,    -0.75);
    dJointSetHinge2Param (joint[0], dParamHiStop,     0.75);
    dJointSetHinge2Param (joint[0], dParamFudgeFactor, 0.1);

    dSpaceCollide (space,0,&nearCallback);
    dWorldStep (world,0.05);

    // remove all contact joints
    dJointGroupEmpty (contactgroup);

  }

  dsSetColor (0,1,1);
  dsSetTexture (DS_WOOD);
  dReal sides[3] = {LENGTH,WIDTH,HEIGHT};
  dsDrawBox (dBodyGetPosition(body[0]),dBodyGetRotation(body[0]),sides);
  dsSetColor (1,1,1);
  for (i=1; i<=3; i++) dsDrawCylinder (dBodyGetPosition(body[i]),
				       dBodyGetRotation(body[i]),0.02f,RADIUS);

  dVector3 ss;
  dGeomBoxGetLengths (ground_box,ss);
  dsDrawBox (dGeomGetPosition(ground_box),dGeomGetRotation(ground_box),ss);

}


int main (int argc, char **argv)
{
  int i;
  dMass m;

  // setup pointers to drawstuff callback functions
  dsFunctions fn;
  fn.version = DS_VERSION;
  fn.start = &start;
  fn.step = &simLoop;
  fn.command = &command;
  fn.stop = 0;
  fn.path_to_textures = DRAWSTUFF_TEXTURE_PATH;

  // create world
  dInitODE2(0);
  world = dWorldCreate();
  space = dHashSpaceCreate (0);
  contactgroup = dJointGroupCreate (0);
  dWorldSetGravity (world,0,0,-0.5);
  ground = dCreatePlane (space,0,0,1,0);

  // chassis body
  body[0] = dBodyCreate (world);
  dBodySetPosition (body[0],0,0,STARTZ);
  dMassSetBox (&m,1,LENGTH,WIDTH,HEIGHT);
  dMassAdjust (&m,CMASS);
  dBodySetMass (body[0],&m);
  box[0] = dCreateBox (0,LENGTH,WIDTH,HEIGHT);
  dGeomSetBody (box[0],body[0]);

  // wheel bodies
  for (i=1; i<=3; i++) {
    body[i] = dBodyCreate (world);
    dQuaternion q;
    dQFromAxisAndAngle (q,1,0,0,M_PI*0.5);
    dBodySetQuaternion (body[i],q);
    dMassSetSphere (&m,1,RADIUS);
    dMassAdjust (&m,WMASS);
    dBodySetMass (body[i],&m);
    sphere[i-1] = dCreateSphere (0,RADIUS);
    dGeomSetBody (sphere[i-1],body[i]);
  }
  dBodySetPosition (body[1],0.5*LENGTH,0,STARTZ-HEIGHT*0.5);
  dBodySetPosition (body[2],-0.5*LENGTH, WIDTH*0.5,STARTZ-HEIGHT*0.5);
  dBodySetPosition (body[3],-0.5*LENGTH,-WIDTH*0.5,STARTZ-HEIGHT*0.5);

  // front and back wheel hinges
  for (i=0; i<3; i++) {
    joint[i] = dJointCreateHinge2 (world,0);
    dJointAttach (joint[i],body[0],body[i+1]);
    const dReal *a = dBodyGetPosition (body[i+1]);
    dJointSetHinge2Anchor (joint[i],a[0],a[1],a[2]);
    dJointSetHinge2Axes (joint[i], zunit, yunit);
  }

  // set joint suspension
  for (i=0; i<3; i++) {
    dJointSetHinge2Param (joint[i],dParamSuspensionERP,0.4);
    dJointSetHinge2Param (joint[i],dParamSuspensionCFM,0.8);
  }

  // lock back wheels along the steering axis
  for (i=1; i<3; i++) {
    // set stops to make sure wheels always stay in alignment
    dJointSetHinge2Param (joint[i],dParamLoStop,0);
    dJointSetHinge2Param (joint[i],dParamHiStop,0);
    // the following alternative method is no good as the wheels may get out
    // of alignment:
    //   dJointSetHinge2Param (joint[i],dParamVel,0);
    //   dJointSetHinge2Param (joint[i],dParamFMax,dInfinity);
  }

  // create car space and add it to the top level space
  car_space = dSimpleSpaceCreate (space);
  dSpaceSetCleanup (car_space,0);
  dSpaceAdd (car_space,box[0]);
  dSpaceAdd (car_space,sphere[0]);
  dSpaceAdd (car_space,sphere[1]);
  dSpaceAdd (car_space,sphere[2]);

  // environment
  ground_box = dCreateBox (space,2,1.5,1);
  dMatrix3 R;
  dRFromAxisAndAngle (R,0,1,0,-0.15);
  dGeomSetPosition (ground_box,2,0,-0.34);
  dGeomSetRotation (ground_box,R);

  // run simulation
  dsSimulationLoop (argc,argv,352,288,&fn);

  dGeomDestroy (box[0]);
  dGeomDestroy (sphere[0]);
  dGeomDestroy (sphere[1]);
  dGeomDestroy (sphere[2]);
  dJointGroupDestroy (contactgroup);
  dSpaceDestroy (space);
  dWorldDestroy (world);
  dCloseODE();
  return 0;
}