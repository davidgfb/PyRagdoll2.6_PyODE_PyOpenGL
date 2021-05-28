/*
buggy with suspension.
this also shows you how to use geom groups.
*/

#include <ode/ode.h>
#include <drawstuff/drawstuff.h>
#include "Semaforo.cpp"

float LENGTH = 5.0, 
	WIDTH = 2.0, 
	HEIGHT = 1.0, 
	RADIUS = 0.5, 
	STARTZ = 2.0, //altura desde la que empieza
	CMASS = 10, 
	WMASS = 10.0,
	v = 10.0,
	fMax = 3.0; //aceleracion
dWorldID world;
dSpaceID space;
dBodyID body[4];
dJointID joint[3];	// joint[0] is the front wheel
dJointGroupID contactgroup;
dGeomID ground,
	box[1],
	sphere[3],
	ground_box;
dSpaceID car_space;
dReal speed=0,
	steer=0;	
int traccion = 3; //0 sin traccion, 1 traccion delantera, 2 traccion trasera, 3 traccion cuatro ruedas
dVector3 yunit = {0, 1, 0}, 
	zunit = {0, 0, 1};

//array de dGeomsIDs
void nearCallback(void *data, dGeomID o1, dGeomID o2) {
  int g1 = (o1 == ground || o1 == ground_box),
	  g2 = (o2 == ground || o2 == ground_box),
	  N = 10;
	dContact contact[N];
	  	  
  if (g1 ^ g2) {
	  int n = dCollide(o1,o2,N,&contact[0].geom,sizeof(dContact));
	  
	  if (n > 0) {
		for (int pos = 0; pos < n; pos++) {
		  contact[pos].surface.mode = dContactSlip1 + dContactSlip2 + dContactSoftERP + dContactSoftCFM + dContactApprox1;
		  contact[pos].surface.mu = dInfinity;
		  contact[pos].surface.slip1 = 0.1;
		  contact[pos].surface.slip2 = 0.1;
		  contact[pos].surface.soft_erp = 0.5;
		  contact[pos].surface.soft_cfm = 0.3;
		  dJointID c = dJointCreateContact(world, contactgroup, &contact[pos]);
		  dJointAttach (c, dGeomGetBody(contact[pos].geom.g1),
						   dGeomGetBody(contact[pos].geom.g2));
		}
	  }
  }
}

void start() {
  dAllocateODEDataForThread(dAllocateMaskAll);

  float pos[3] = {5.0, -5.0, 5.0},
  	rot[3] = {121.0, -45.0, 0.0};
  dsSetViewpoint (pos, rot);
}

void command (int cmd) {
  switch (cmd) {
	  case 'w': 
	  case 'W':
		speed = v;
		break;
	  case 's': 
	  case 'S':
		speed = -v;
		break;
	  	  
	  case 'a':
	  case 'A':
		steer -= 0.5;
		break;
	  case 'd':
	  case 'D':
		steer += 0.5;
		break;
	  
	  case ' ':
		speed = 0;
		steer = 0;
		break;
	}
}



Semaforo semaforo;

void tracciona() {
  	for (int pos = 0; pos < traccion; pos++) {
  		if (traccion == 1 || traccion == 2 && pos > 0 || traccion == 3) {
			dJointSetHinge2Param (joint[pos],dParamVel2,-speed); 
			dJointSetHinge2Param (joint[pos],dParamFMax2,fMax);
		}
    }
}

void gira() {
    dReal v = 10.0 * (steer - dJointGetHinge2Angle1(joint[0]));
    
    if (v < -1.0) {
    	v = -1.0;
    }
    
    if (v > 1.0) {
    	v = 1.0;
    }
            
    dJointSetHinge2Param (joint[0],dParamVel,v);
    dJointSetHinge2Param (joint[0],dParamFMax, 0.2);
    dJointSetHinge2Param (joint[0],dParamLoStop, -0.75);
    dJointSetHinge2Param (joint[0],dParamHiStop, 0.75);
    dJointSetHinge2Param (joint[0],dParamFudgeFactor, 0.1);
}

void simLoop (int pause) {
  if (!pause) {  
  	tracciona();

	gira();

    dSpaceCollide (space, 0, &nearCallback);
    dWorldStep (world, 0.1); //camara lenta

    dJointGroupEmpty (contactgroup);
  }

  dsSetColor(0, 1,1);
  dReal sides[3] = {LENGTH,WIDTH,HEIGHT};
   dsDrawBox(dBodyGetPosition(body[0]),
     		 dBodyGetRotation(body[0]),sides);
   dsSetColor (1,1,1);
  	
    semaforo.actualiza();
  
	for (int pos = 1; pos < 4; pos++) {
	  	dsDrawCylinder (dBodyGetPosition(body[pos]),
						dBodyGetRotation(body[pos]),
						0.02, RADIUS);
	}

  dVector3 ss;
  dGeomBoxGetLengths (ground_box,ss);
  dsDrawBox(dGeomGetPosition(ground_box),
  			dGeomGetRotation(ground_box),ss);
}

int main (int argc, char **argv) {
  dMass m;

  dsFunctions fn;
  fn.version = DS_VERSION;
  fn.start = &start;
  fn.step = &simLoop;
  fn.command = &command;
  fn.stop = 0;
  fn.path_to_textures = DRAWSTUFF_TEXTURE_PATH;

  dInitODE2(0);
  world = dWorldCreate();
  space = dHashSpaceCreate(0);
  contactgroup = dJointGroupCreate (0);
  dWorldSetGravity (world, 0,0,-0.5); //
  ground = dCreatePlane (space, 0,0,1,0);

  body[0] = dBodyCreate (world);
  dBodySetPosition (body[0], 0,0, STARTZ);
  dMassSetBox (&m, 1, LENGTH,WIDTH,HEIGHT);
  dMassAdjust (&m,CMASS);
  dBodySetMass (body[0],&m);
  box[0] = dCreateBox (0, LENGTH,WIDTH,HEIGHT);
  dGeomSetBody (box[0],body[0]);

  for (int pos = 1; pos < 4; pos++) {
    body[pos] = dBodyCreate (world);
    dQuaternion q;
    dQFromAxisAndAngle (q, 1,0,0, M_PI / 2.0);
    dBodySetQuaternion (body[pos],q);
    dMassSetSphere (&m, 1, RADIUS);
    dMassAdjust (&m,WMASS);
    dBodySetMass (body[pos],&m);
    sphere[pos - 1] = dCreateSphere (0, RADIUS);
    dGeomSetBody (sphere[pos - 1],body[pos]);
  }
  
  dBodySetPosition (body[1], LENGTH / 2.0, 0, STARTZ - HEIGHT / 2.0);
  dBodySetPosition (body[2], LENGTH / -2.0, WIDTH / 2.0, STARTZ - HEIGHT / 2.0);
  dBodySetPosition (body[3], LENGTH / -2.0, WIDTH / -2.0, STARTZ - HEIGHT / 2.0);

  for (int pos = 0; pos < 3; pos++) {
    joint[pos] = dJointCreateHinge2 (world, 0);
    dJointAttach (joint[pos],body[0],body[pos + 1]);
    const dReal *a = dBodyGetPosition (body[pos + 1]);
    dJointSetHinge2Anchor (joint[pos],a[0],a[1],a[2]);
    
    dJointSetHinge2Axes   (joint[pos], zunit, yunit);
  }

  for (int pos = 0; pos < 3; pos++) {
    dJointSetHinge2Param (joint[pos],dParamSuspensionERP, 0.4);
    dJointSetHinge2Param (joint[pos],dParamSuspensionCFM, 0.8);
  }

  for (int pos = 1; pos < 3; pos++) {
    dJointSetHinge2Param (joint[pos],dParamLoStop, 0);
    dJointSetHinge2Param (joint[pos],dParamHiStop, 0);

  }

  car_space = dSimpleSpaceCreate (space);
  dSpaceSetCleanup (car_space, 0);
  dSpaceAdd (car_space,box[0]);
  
  for (int pos = 0; pos < 3; pos++) {
	  dSpaceAdd (car_space,sphere[pos]);
	}

  ground_box = dCreateBox(space, 5,5,5);
  dMatrix3 R;
  dRFromAxisAndAngle (R, 0,1,0,-0.15);
  dGeomSetPosition (ground_box, 2,0,-2.0);
  dGeomSetRotation (ground_box,R);

  dsSimulationLoop (argc,argv, 800,600, &fn);

  return 0;
}




