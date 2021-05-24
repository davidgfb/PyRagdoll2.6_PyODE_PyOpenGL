#include <ode/ode.h>
#include <drawstuff/drawstuff.h>
#include "texturepath.h"

struct MyObject {
	dBodyID body;			
	dGeomID geom[1];		
};

int show_contacts = 0;	

dWorldID world;
dSpaceID space;
MyObject obj[1];
dJointGroupID contactgroup;

dReal dVector3R[3];

static void nearCallback(void*, dGeomID o1, dGeomID o2) {
	dBodyID b1 = dGeomGetBody( o1 ),
	b2 = dGeomGetBody( o2 );
	
	if (not( b1 && b2 &&  dAreConnectedExcluding( b1,b2,dJointTypeContact ))) {

		dContact contact[1];   
		
		dContact contacto = contact[0];
		
			contacto.surface.mode = dContactBounce + dContactSoftCFM; 
			contacto.surface.mu = dInfinity;
			contacto.surface.mu2 = 0;
			contacto.surface.bounce = 0.1;
			contacto.surface.bounce_vel = 0.1;
			contacto.surface.soft_cfm = 0.01;
		
		dCollide(o1, o2, 1, &contact[0].geom, sizeof(dContact));
		//&contact[0] != &contacto
		dMatrix3 RI;
		dRSetIdentity(RI);
		dReal ss[3] = {0.02, 0.02, 0.02};
			
		dJointAttach(dJointCreateContact(world, contactgroup, contact), b1, b2);		
	}
}

int i = 0;

static void start() {
	static float xyz[3] = {1.0, 0.0, 1.0},
	hpr[3] = {180.0, -45.0, 0.0};
	dsSetViewpoint(xyz, hpr);
		
	dReal sides[3];
	dMass m;
		
	obj[i].body = dBodyCreate( world );
		
	for (int k = 0; k < 3; k++) {
		sides[k] = 0.2;
	}

	dMatrix3 R;
	
	dBodySetPosition(obj[i].body, 0, 0, 1);			
	dRFromAxisAndAngle(R, 0.0, 0.0, 0.0, 0.0 );				
	dBodySetRotation(obj[i].body, R);
	dBodySetData(obj[i].body, (void*)(dsizeint) i);

	sides[0] /= 2.0;
	dMassSetCapsule(&m, 5.0, 3, sides[0], sides[1]);
	obj[i].geom[0] =  dCreateCapsule(space, sides[0], sides[1]);
      
	dGeomSetBody(obj[i].geom[0], obj[i].body );
	dBodySetMass(obj[i].body, &m);
}

void drawGeom( dGeomID g, const dReal *pos, const dReal *R, int show_aabb ) {
		dReal radius,length;
		dGeomCapsuleGetParams( g,&radius,&length );
		dsDrawCapsule(dGeomGetPosition(g), dGeomGetRotation(g), length, radius);		
}

static void simLoop(int pause) {
	dSpaceCollide(space, 0, &nearCallback);
	dWorldQuickStep(world, 0.05);

	for (int j = 0; j < 2; j++) {
		dSpaceGetGeom(space, j);
	}

	dJointGroupEmpty(contactgroup);
	drawGeom(obj[0].geom[0], 0, 0, 0);
}


int main(int argc, char **argv){
	dsFunctions fn;
	
	fn.version = DS_VERSION;
	fn.start = &start;
	fn.step = &simLoop;
	fn.stop = 0;
	fn.path_to_textures = DRAWSTUFF_TEXTURE_PATH;

	dInitODE2(0);
	world = dWorldCreate();

	space = dSimpleSpaceCreate(0);
	contactgroup = dJointGroupCreate(0);
	dWorldSetGravity(world, 0, 0, -0.5 );
	dWorldSetCFM(world, 1e-5);
	dCreatePlane(space, 0, 0, 1,0);

	dsSimulationLoop(argc, argv, 1000, 1000, &fn);

	dJointGroupDestroy(contactgroup);
	dSpaceDestroy(space);
	dWorldDestroy(world);
	dCloseODE();
	
	return 0;
}

 	  	 

