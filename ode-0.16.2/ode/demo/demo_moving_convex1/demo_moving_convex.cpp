#include <ode/ode.h>
#include <drawstuff/drawstuff.h>
#include "texturepath.h"
#include "bunny_geom.h"
#include "convex_bunny_geom.h"

#ifdef _MSC_VER
	#pragma warning(disable:4244 4305)  // for VC++, no precision loss complaints
#endif

// select correct drawing functions
#ifdef dDOUBLE
	#define dsDrawBox dsDrawBoxD
	#define dsDrawSphere dsDrawSphereD
	#define dsDrawCylinder dsDrawCylinderD
	#define dsDrawCapsule dsDrawCapsuleD
	#define dsDrawLine dsDrawLineD
	#define dsDrawTriangle dsDrawTriangleD
	#define dsDrawConvex dsDrawConvexD
#endif

// some constants
#define NUM 200			// max number of objects
#define DENSITY (5.0)		// density of all objects
#define GPB 3			// maximum number of geometries per body
#define MAX_CONTACTS 64		// maximum number of contact points per body

/*
int NUM = 200,
GPB = 3,
MAX_CONTACTS = 64;
float DENSITY = 5.0;
*/

// dynamics and collision objects
struct MyObject {
	dBodyID body;			// the body
	dGeomID geom[GPB];		// geometries representing this body
};

static int num = 0,		// number of objects in simulation
nextobj = 0,		// next object to recycle if num==NUM
selected = -1,	// selected object
show_aabb = 0,	// show geom AABBs?
show_contacts = 0,	// show contact points?
random_pos = 1;	// drop objects from random position?

static dWorldID world;
static dSpaceID space;
static MyObject obj[NUM];
static dJointGroupID contactgroup;

typedef dReal dVector3R[3];
// this is called by dSpaceCollide when two objects in space are
// potentially colliding.

static void nearCallback( void *, dGeomID o1, dGeomID o2 ){
	// if (o1->body && o2->body) return;

	// exit without doing anything if the two bodies are connected by a joint
	dBodyID b1 = dGeomGetBody( o1 ),
	b2 = dGeomGetBody( o2 );
	
	if (not( b1 && b2 &&  dAreConnectedExcluding( b1,b2,dJointTypeContact ))) {

		dContact contact[MAX_CONTACTS];   // up to MAX_CONTACTS contacts per box-box
		
		for (int i = 0; i < MAX_CONTACTS; i++){
			contact[i].surface.mode = dContactBounce | dContactSoftCFM;
			contact[i].surface.mu = dInfinity;
			contact[i].surface.mu2 = 0;
			contact[i].surface.bounce = 0.1;
			contact[i].surface.bounce_vel = 0.1;
			contact[i].surface.soft_cfm = 0.01;
		}
		
		if (int numc =  dCollide(o1, o2, MAX_CONTACTS, &contact[0].geom, sizeof(dContact))) {
			dMatrix3 RI;
			dRSetIdentity(RI);
			const dReal ss[3] = {0.02, 0.02, 0.02};
			
			for (int i = 0; i < numc; i++){
				dJointID c =  dJointCreateContact(world,contactgroup,contact+i );
				dJointAttach(c, b1, b2);
				
				if (show_contacts) { 
					dsDrawBox(contact[i].geom.pos, RI, ss);
				}
			}
		}
	}
}

/*
char locase( char c ) { //interesante
	char cMinusculas = c;
	
	if (c >= 'A' && c <= 'Z') {
		cMinusculas = c - ( 'a'-'A' );
	} 
	
	return cMinusculas;
}
*/

// start simulation - set viewpoint
static void start(){
	dAllocateODEDataForThread(dAllocateMaskAll);

	static float xyz[3] = {2.0, -1.0, 2.0},
	hpr[3] = {126.0, -17.0, 0.0};
	dsSetViewpoint(xyz, hpr);
	
	int i = 0;
	dReal sides[3];
	dMass m;
	
	if (num < NUM) {
		i = num;
		num++;
		
	} else {
		i = nextobj;
		nextobj++;
			
		if (nextobj >= num) {
			nextobj = 0;
		}

		// destroy the body and geoms for slot i
		dBodyDestroy(obj[i].body);
			
		for (int k = 0; k < GPB; k++){
			if (obj[i].geom[k]) {
					dGeomDestroy(obj[i].geom[k]);
			}
		}
			
		memset(&obj[i], 0, sizeof(obj[i])); //
	}

	obj[i].body = dBodyCreate( world );
		
	for (int k = 0; k < 3; k++) {
		sides[k] = dRandReal() / 2.0 + 0.1;
	}

	dMatrix3 R;
		
	if (random_pos) {
			dBodySetPosition(obj[i].body, dRandReal() * 2 - 1, dRandReal() * 2 - 1, dRandReal() + 3);
			dRFromAxisAndAngle(R, dRandReal() * 2.0 - 1.0, dRandReal() * 2.0 - 1.0, dRandReal() * 2.0 - 1.0, dRandReal() * 10.0 - 5.0 );
		
		} else {
			dReal maxheight = 0;
			
			for (int k = 0; k < num; k++ ){
				const dReal *pos =  dBodyGetPosition(obj[k].body);
				
				if (pos[2] > maxheight) {
					maxheight = pos[2];
				}
			}
			
			dBodySetPosition(obj[i].body, 0, 0, maxheight + 1);
			dRFromAxisAndAngle(R, 0, 0, 1, dRandReal() * 10.0 - 5.0);
		}
		
		dBodySetRotation(obj[i].body, R);
		dBodySetData(obj[i].body, (void*)(dsizeint) i);

		sides[0] /= 2.0;
		dMassSetCapsule(&m, DENSITY, 3, sides[0], sides[1]);
		obj[i].geom[0] =  dCreateCapsule(space, sides[0], sides[1]);
      
		for (int k = 0; k < GPB; k++) {
			if (obj[i].geom[k]) {
			 	dGeomSetBody(obj[i].geom[k], obj[i].body );
			}
		}

		dBodySetMass(obj[i].body, &m);
}

// draw a geom
void drawGeom( dGeomID g, const dReal *pos, const dReal *R, int show_aabb ) {
	if (g) {		
		if (!pos) {
			pos = dGeomGetPosition(g);
		}
		
		if (!R) {
			R = dGeomGetRotation(g);
		}

		int type = dGeomGetClass(g);
		
		if (type == dCapsuleClass ) {
			dReal radius,length;
			dGeomCapsuleGetParams( g,&radius,&length );
			dsDrawCapsule(pos,R,length,radius);		
		} 
		
		if ( show_aabb ) {
			// draw the bounding box for this geom
			dReal aabb[6];
			dGeomGetAABB( g,aabb );
			dVector3 bbpos;
			
			for (int i = 0; i < 3; i++) {
				bbpos[i] = ( aabb[i * 2] + aabb[i * 2 + 1] ) / 2.0;
			}
			
			dVector3 bbsides;
			
			for (int j = 0; j < 3; j++) {
				bbsides[j] = aabb[j * 2 + 1] - aabb[j * 2];
			}
			
			dMatrix3 RI;
			dRSetIdentity( RI );
			dsSetColorAlpha(1, 0, 0, 0.5);
			dsDrawBox( bbpos,RI,bbsides );
		}
	}
}

// simulation loop
static void simLoop(int pause) {
	dsSetColor(0, 0, 2);
	dSpaceCollide(space, 0, &nearCallback);

	if (!pause) {
		dWorldQuickStep(world, 0.05);
	}

	for (int j = 0; j < dSpaceGetNumGeoms(space); j++) {
		dSpaceGetGeom(space, j);
	}

	// remove all contact joints
	dJointGroupEmpty(contactgroup);

	dsSetColor(1, 1, 0);
	dsSetTexture(DS_WOOD);
	for (int i = 0; i < num; i++){
		for (int j = 0; j < GPB; j++){
			if (obj[i].geom[j]){
				if (i == selected){
					dsSetColor(0, 0.7, 1);
				
				} else if (!dBodyIsEnabled(obj[i].body)) {
					dsSetColor(1,0,0);
				
				} else {
					dsSetColor(1,1,0);
				}

				drawGeom(obj[i].geom[j], 0, 0, show_aabb);
			}
		}
	}
}


int main(int argc, char **argv){
	// setup pointers to drawstuff callback functions
	dsFunctions fn;
	fn.version = DS_VERSION;
	fn.start = &start;
	fn.step = &simLoop;
	//fn.command = &command;
	fn.stop = 0;
	fn.path_to_textures = DRAWSTUFF_TEXTURE_PATH;

	// create world
	dInitODE2(0);
	world = dWorldCreate();

	space = dSimpleSpaceCreate(0);
	contactgroup = dJointGroupCreate(0);
	dWorldSetGravity(world, 0, 0, -0.5 );
	dWorldSetCFM(world, 1e-5);
	dCreatePlane(space, 0, 0, 1,0);
	memset(obj, 0, sizeof(obj)); //look

	// run simulation
	dsSimulationLoop(argc, argv, 1000, 1000, &fn);

	dJointGroupDestroy(contactgroup);
	dSpaceDestroy(space);
	dWorldDestroy(world);
	dCloseODE();
	
	return 0;
}

 	  	 
