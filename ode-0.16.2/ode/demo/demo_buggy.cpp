/*
buggy with suspension.
this also shows you how to use geom groups.
*/
#include <ode/ode.h>
#include <drawstuff/drawstuff.h>
#include "texturepath.h"

//////////// variables ////////////////////
float   LENGTH = 0.7, 
	WIDTH = 0.5, 
	HEIGHT = 0.2, 
	RADIUS = 0.18, 
	STARTZ = 0.5, 
	CMASS = 1.0, 
	WMASS = 0.2; // some constants

static const dVector3   yunit = {0, 1, 0}, 
			zunit = {0, 0, 1};

static dWorldID world; // dynamics and collision objects (chassis, 3 wheels, environment)
static dSpaceID space, car_space;
static dBodyID body[4];
static dJointID joint[3]; // array3d de dJointIDs
static dJointGroupID contactgroup;
static dGeomID ground, box[1], sphere[3], ground_box;
static dReal    speed = 0, 
		steer = 0; // things that the user controls // user commands
///////////////////////////////////////////

//////////////// variables ragdoll /////////////
//rotation directions are named by the third (z-axis) row of the 3x3 matrix, because ODE capsules are oriented along the z-axis
float rightRot[9] = {0.0, 0.0, -1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0},
      leftRot[9]  = {0.0, 0.0, 1.0, 0.0, 1.0, 0.0, -1.0, 0.0, 0.0},
      upRot[9]    = {1.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 1.0, 0.0},
      downRot[9]  = {1.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 1.0, 0.0},
      bkwdRot[9]  = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0},

//axes used to determine constrained joint rotations
rightAxis[3] = { 1.0,  0.0,  0.0},
leftAxis[3]  = {-1.0,  0.0,  0.0},
upAxis[3]    = { 0.0,  1.0,  0.0},
downAxis[3]  = { 0.0, -1.0,  0.0},
bkwdAxis[3]  = { 0.0,  0.0,  1.0},
fwdAxis[3]   = { 0.0,  0.0, -1.0},

UPPER_ARM_LEN =  0.3,
FORE_ARM_LEN  = 0.25,
HAND_LEN      = 0.13, // wrist to mid-fingers only
FOOT_LEN      = 0.18, //ankles to base of ball of foot only
HEEL_LEN      = 0.05,

BROW_H     = 1.68,
MOUTH_H    = 1.53,
NECK_H     =  1.5,
SHOULDER_H = 1.37,
CHEST_H    = 1.35,
HIP_H      = 0.86,
KNEE_H     = 0.48,
ANKLE_H    = 0.08,

SHOULDER_W = 0.41,
CHEST_W    = 0.36, // actually wider, but we want narrower than shoulders (esp. with large radius)
LEG_W      = 0.28, // between middles of upper legs
PELVIS_W   = 0.25, // actually wider, but we want smaller than hip width;

R_SHOULDER_POS[3] = {(float) (SHOULDER_W / -2.0), SHOULDER_H, 0.0},
L_SHOULDER_POS[3] = {(float) (SHOULDER_W / 2.0), SHOULDER_H, 0.0},

R_ELBOW_POS[3],
L_ELBOW_POS[3],
R_WRIST_POS[3],
L_WRIST_POS[3],
R_FINGERS_POS[3],
L_FINGERS_POS[3],

R_HIP_POS[3]   = {(float) (LEG_W / -2.0),   HIP_H, 0.0},
L_HIP_POS[3]   = {(float) (LEG_W /  2.0),   HIP_H, 0.0},
R_KNEE_POS[3]  = {(float) (LEG_W / -2.0),  KNEE_H, 0.0},
L_KNEE_POS[3]  = {(float) (LEG_W /  2.0),  KNEE_H, 0.0},
R_ANKLE_POS[3] = {(float) (LEG_W / -2.0), ANKLE_H, 0.0},
L_ANKLE_POS[3] = {(float) (LEG_W /  2.0), ANKLE_H, 0.0},

R_HEEL_POS[3],
L_HEEL_POS[3],
R_TOES_POS[3],
L_TOES_POS[3];


////////////////////////////////////////////////

///////////// ragdoll //////////////////
float sign(float x) { //sentido bool
	//Returns 1.0 if x is positive, -1.0 if x is negative or zero.
	float valor = -1.0; // (x <= 0.0), x E (-inf, 0]
	
	if (x > 0.0) { // (0, inf)
		valor = 1.0;
	} 
	
	return valor;
}

float len3(float v[3]) { 
	//Returns the length of 3-vector v.
	float x = v[0], y = v[1], z = v[2];
	
	return sqrt(x * x + y * y + z * z);
}

//error sizeof(puntero) // tamaño
void imprimeArray(float* array, int t) { //float star = tipo array	
	for (int posicion = 0; posicion < t; posicion++) { 
		printf("%f ", array[posicion]); //array);
	}	
}

void neg3(float v[3]) { //no devuelve nada haz una copia antes
	//the negation of 3-vector v.
	for (int posicion = 0; posicion < 3; posicion++) {
		v[posicion] *= -1.0; 
	}
}

void add3(float a[3], float b[3]) { //modifica a haz una copia
	//the sum of 3-vectors a and b.
	a[0] += b[0];
	a[1] += b[1];
	a[2] += b[2];
}

void sub3(float a[3], float b[3]) { //modifica a haz una copia
	//the difference between 3-vectors a and b.
	a[0] -= b[0];
	a[1] -= b[1];
	a[2] -= b[2];
}

void mul3(float v[3], float s) { //haz una copia
	//3-vector v multiplied by scalar s.
	v[0] *= s;
	v[1] *= s;
	v[2] *= s;
}

void div3(float v[3], float s) {
	//3-vector v divided by scalar s.
	v[0] /= s;
	v[1] /= s;
	v[2] /= s;
}


float dist3(float a[3], float b[3]) {
	//Returns the distance between point 3-vectors a and b.		
	sub3(a, b);
		
	return len3(a);
}

void norm3(float v[3]) {
	//the unit length 3-vector parallel to 3-vector v.		
	float l = len3(v); //modifica v
	
	if (l > 0.0) { //funcion bool
		div3(v, l);
	} else {
		v[0] = 0.0;
		v[1] = 0.0;
		v[2] = 0.0;
	}
}

float dot3(float a[3], float b[3]) {
	//Returns the dot product of 3-vectors a and b.
	float xa = a[0], ya = a[1], za = a[2], 
	      xb = b[0], yb = b[1], zb = b[2];
	      
	return xa * xb + ya * yb + za * zb;
}

void cross(float a[3], float b[3]) { //por que no va?
	//the cross product of 3-vectors a and b.	
	float xa = a[0], ya = a[1], za = a[2], 
	      xb = b[0], yb = b[1], zb = b[2];
	      
	a[0] = ya * zb - yb * za;
	a[1] = za * xb - zb * xa;
	a[2] = xa * yb - xb * ya;
}

void project3(float v[3], float d[3]) {
	//projection of 3-vector v onto unit 3-vector d.
	float v1[3] = {v[0], v[1], v[2]}; //copia de v 
	
	norm3(v1); //modifica v1
		
	mul3(v, dot3(v1, d));
}

float acosdot3(float a[3], float b[3]) {
	//Returns the angle between unit 3-vectors a and b.
	float x = dot3(a, b), valor = 0.0; // valor E (0, inf)
	
	if (x < -1.0) { //llama
		valor = M_PI;
	} else {
		if (x == 0.0) {
			valor = acos(x); //acos(0)
		}
	}
	
	return valor;
}

void invert3x3(float m[9]) {
	//the inversion (transpose) of 3x3 rotation matrix m.
	float            m1 = m[1], m2 = m[2], //m0
	      m3 = m[3], m4 = m[4], m5 = m[5], 
	      m6 = m[6], m7 = m[7], m8 = m[8];
	
	m[1] = m3;
	m[2] = m6;
	
	m[3] = m1;
	m[4] = m4;
	m[5] = m7;
	
	m[6] = m2;
	m[7] = m5;
	m[8] = m8;
}



void rotate3(float m[9], float v[3]) { //guarda en v! no en m!
	//the rotation of 3-vector v by 3x3 (row major) matrix m.
	float xv = v[0], yv = v[1], zv = v[2];
	
	v[0] = xv * m[0] + yv * m[1] + zv * m[2];
	v[1] = xv * m[3] + yv * m[4] + zv * m[5];
	v[2] = xv * m[6] + yv * m[7] + zv * m[8];
	
	//v = {xv * m[0] + yv * m[1] + zv * m[2],
	//     xv * m[3] + yv * m[4] + zv * m[5],
 	//     xv * m[6] + yv * m[7] + zv * m[8]};
}

void zaxis(float m[9]) {
	//the z-axis vector from 3x3 (row major) rotation matrix m.
	m[0] = m[2];
	m[1] = m[5];
	m[2] = m[8];
}

void calcRotMatrix(float axis[9], float angle) {
	//the row-major 3x3 rotation matrix defining a rotation around axis by	angle.
	float cosTheta = cos(angle), 
	      sinTheta = sin(angle),
	      t = 1.0 - cosTheta,
	      
	      axis0 = axis[0],
	      axis1 = axis[1],
	      axis2 = axis[2],
	      
	      axis3 = axis[3],
	      axis4 = axis[4],
	      axis5 = axis[5],
	      
	      axis6 = axis[6],
	      axis7 = axis[7],
	      axis8 = axis[8];
	
	axis[0] = t * axis0 * axis0 + cosTheta;
	axis[1] = t * axis0 * axis1 - sinTheta * axis2;
	axis[2] = t * axis0 * axis2 + sinTheta * axis1;
	
	axis[3] = t * axis0 * axis1 + sinTheta * axis2;
	axis[4] = t * axis1 * axis1 + cosTheta;
	axis[5] = t * axis1 * axis2 - sinTheta * axis0;
	
	axis[6] = t * axis0 * axis2 - sinTheta * axis1;
	axis[7] = t * axis1 * axis2 + sinTheta * axis0;
	axis[8] = t * axis2 * axis2 + cosTheta;
}


void makeOpenGLMatrix(float r[16], float p[3]) { //no puede ser float r[9] -> r[16]
	//an OpenGL compatible (column-major, 4x4 homogeneous) transformation	matrix from ODE compatible (row-major, 3x3) rotation matrix r and position	vector p.
	float            r1 = r[1], r2 = r[2], 
	      r3 = r[3], r4 = r[4], r5 = r[5], 
	      r6 = r[6], r7 = r[7], r8 = r[8];
	
	//r = {};
	r[1] = r3;
	r[2] = r6;
	
	r[3] = 0.0;	
	r[4] = r1;
	r[5] = r4;	
	
	r[6] = r7;
	r[7] = 0.0;	
	r[8] = r2; //hasta aqui es r
	
	r[9] = r5; //a partir de aqui es r+ extendido
	r[10] = r8;
	r[11] = 0.0;
	
	r[12] = p[0];
	r[13] = p[1];
	r[14] = p[2];
	r[15] = 1.0;	
}

void asigna_Array(float *array, float *array1, int t) {
	for (int pos = 0; pos < t; pos++) {
		array[pos] = array1[pos];
	}
}

//obsoleto?
void asigna_Array3(float array[3], float array1[3]) {
	asigna_Array(array,array1,3);
}

//////////////////////////////////////

dJointID getJuntaRuedaDelantera() {
	return joint[0];
}

static void nearCallback (void *, dGeomID o1, dGeomID o2) {
	  int i = 0, n = 0;

	  bool g1 = (o1 == ground || o1 == ground_box), 	  // only collide things with the ground
	  g2 = (o2 == ground || o2 == ground_box);
	  
	  if (g1 ^ g2) {
		  const int N = 10;
		  dContact contact[N];
		  n = dCollide (o1, o2, N, &contact[0].geom, sizeof(dContact));
		  
		  int sumaContactos = dContactSlip1   + dContactSlip2   +
				      dContactSoftERP + dContactSoftCFM + dContactApprox1;
		  
			    for (i = 0; n>0 && i < n; i++) {
			      contact[i].surface.mode = sumaContactos;
			      
			      contact[i].surface.mu = dInfinity;
			      
			      contact[i].surface.slip1    = 0.1;
			      contact[i].surface.slip2    = 0.1;
			      contact[i].surface.soft_erp = 0.5;
			      contact[i].surface.soft_cfm = 0.3;
			      
			      dJointID c = dJointCreateContact (world, contactgroup, &contact[i]);
			      dJointAttach (c, dGeomGetBody(contact[i].geom.g1),
					       dGeomGetBody(contact[i].geom.g2));
			    }

	  }
}

static void start() { // start simulation - set viewpoint
	  dAllocateODEDataForThread(dAllocateMaskAll);

	  static float pos[3] = {    1,      -1,   1}, // array3d de floats 
	               rot[3] = {121.0,   -27.5, 0.0};
	  
	  //camara
	  dsSetViewpoint (pos, rot);
	  printf ("Press:\t'a' to increase speed.\n"
		  	"\t'z' to decrease speed.\n"
			"\t',' to steer left.\n"
			"\t'.' to steer right.\n"
			"\t' ' to reset speed and steering.\n");
}

static void command (int cmd) {  // lower 
					
}

int ciclosEspera = 20, cicloActual = 0;

////////////////// simLoop ////////////////////////
static void simLoop (int pause) { // simulation loop
   	int i = 0;
	float vMax = 0.1;
		  
    	if (!pause) {
         	dJointSetHinge2Param (getJuntaRuedaDelantera(),  dParamVel2, -speed); // motor rueda/s delantera
         	dJointSetHinge2Param (getJuntaRuedaDelantera(), dParamFMax2,    0.1);

    	 	dReal v = steer - dJointGetHinge2Angle1 (getJuntaRuedaDelantera()); // steering
    

    	 	
    	 	dJointSetHinge2Param (getJuntaRuedaDelantera(),         dParamVel,     v);
    	 	dJointSetHinge2Param (getJuntaRuedaDelantera(),        dParamFMax,   0.2);
    	 	dJointSetHinge2Param (getJuntaRuedaDelantera(),      dParamLoStop, -0.75);
    	 	dJointSetHinge2Param (getJuntaRuedaDelantera(),      dParamHiStop,  0.75);
    	 	dJointSetHinge2Param (getJuntaRuedaDelantera(), dParamFudgeFactor,   0.1);

    	 	dSpaceCollide (space,    0, &nearCallback);
    	 	dWorldStep    (world, 0.05);

    	 	dJointGroupEmpty (contactgroup); // remove all contact joints
    	}

  	dsSetColor (0, 1, 1);
  	dsSetTexture (DS_WOOD);
  	dReal sides[3] = {LENGTH, WIDTH, HEIGHT};
  	dsDrawBox (dBodyGetPosition(body[0]), dBodyGetRotation(body[0]), sides);
  	dsSetColor (1, 1, 1);
	
	for (i = 1; i < 4; i++) {
  		 dsDrawCylinder (dBodyGetPosition(body[i]),
				 dBodyGetRotation(body[i]), 0.02, RADIUS);
	}

  	dVector3 ss;
  	dGeomBoxGetLengths (ground_box, ss);
  	dsDrawBox (dGeomGetPosition(ground_box), dGeomGetRotation(ground_box), ss);
}
////////////////// fin simLoop ////////////////////

///////////////////// ragdoll /////////////////////


class Ragdoll {
	//Creates a ragdoll of standard size at the given offset.
	//constructor
	//public: 	
	
	dWorldID world; 
	dSpaceID space;
	float densidad = 0.0,

	      totalMass = 0.0,
		  offset[3] = {0.0, 0.0, 0.0}; //new float[3]
	
	
	Ragdoll(dWorldID w, dSpaceID s, float d, float o[3]) { //no densidad		
		world = w;
		space = s;
		densidad = d;
		asigna_Array3(offset, o); //offset = {o[0],o[1],o[2]};
		
	}
	
	//sobrecarga para tener parametro offset definido por defecto
	Ragdoll(dWorldID world, dSpaceID space, float density) {
		float offset[3] = {0.0, 0.0, 0.0};
		Ragdoll(world,space,density,offset);
	}
	
	void addBody(float p1[3], float p2[3], float radius) {
		//Adds a capsule body between joint positions p1 and p2 and with given radius to the ragdoll.
		add3(p1, offset);
		add3(p2, offset);
		
		//cylinder length not including endcaps, make capsules overlap by half radius at joints
		float cyllen = dist3(p1, p2) - radius; //largo cilindro
	
		dBodyID body = dBodyCreate(world);
	
		dMass m; //?
		
		dMassSetCylinder(&m,densidad,3,radius,cyllen); //dCreateCapsule? //3=eje z //m.
	}
};


class Caja {
	dGeomID ID;

	public:
		Caja(dSpaceID espacio, dReal dx, dReal dy, dReal dz, dMatrix3 R) {
			ID = dCreateBox(space, dx, dy, dz);   
			dGeomSetPosition(ID, 2, 0, -0.34);
  			dGeomSetRotation(ID, R);
		}
		
		dGeomID getID() {
			return ID;
		}
};

class Espacio {
	dSpaceID ID;

	public:
		Espacio(dSpaceID espacio) {
			ID = dSimpleSpaceCreate (space);   	// create car space and add it to the top level space
		}
		
		dSpaceID getID() {
			return ID;
		}
		
		void annade(dSpaceID ID, dGeomID ID_Geom) {
			dSpaceAdd (ID, ID_Geom);
		}
};

class Cuerpo {
	//dReal *posicion;
	dBodyID ID;

	public:
		Cuerpo() {
		
		}
	
		Cuerpo(dWorldID ID_Mundo) {
			ID = dBodyCreate(ID_Mundo);
		}
		
		void setPosicion(dBodyID pID, dReal x, dReal y, dReal z) {					
			dBodySetPosition(pID,x,y,z);
		}
		
		void setMasa(dBodyID ID, dMass *mass) { 
			dBodySetMass(ID,mass);
		}
		
		void setQuaternion(dBodyID ID, dQuaternion q) {
			dBodySetQuaternion(ID,q);
		}
				
		dBodyID getID() {
			return ID;
		}
		
		/*
		dReal* getPosicion() {
			return dBodyGetPosition(ID);
		}
		*/		
};

class Junta {
	public:
		Junta() {
		
		}
		
		void setAnclaBisagra2(dJointID ID, dReal x, dReal y, dReal z) { 
			dJointSetHinge2Anchor (ID,  x,  y, z);
		}
		
		void setEjesBisagra2(dJointID ID, const dReal* eje, const dReal* eje1) {
			dJointSetHinge2Axes(ID,eje,eje1); 
		}
		
		void setParamsBisagra2(dJointID ID, int parameter, dReal value) {
			dJointSetHinge2Param(ID,parameter,value);
		}
};

Cuerpo cuerpos[3];

//clase Coche, Moto...

/////////////// main ///////////////////////
int main (int argc, char **argv) {
	//////////// variables ragdoll ////////////
	float largoAntebrazo[3] = {UPPER_ARM_LEN, 0.0, 0.0};
	sub3(R_SHOULDER_POS, largoAntebrazo);
	asigna_Array3(R_ELBOW_POS, R_SHOULDER_POS);

	add3(L_SHOULDER_POS, largoAntebrazo);	
	asigna_Array3(L_ELBOW_POS, L_SHOULDER_POS);
	
	float largoBrazo[3] = {FORE_ARM_LEN, 0.0, 0.0};
	sub3(R_ELBOW_POS, largoBrazo);
	asigna_Array3(R_WRIST_POS, R_ELBOW_POS);
	
	add3(L_ELBOW_POS, largoBrazo);
	asigna_Array3(L_WRIST_POS, L_ELBOW_POS);
	
	float largoMano[3] = {HAND_LEN, 0.0, 0.0};
	sub3(R_WRIST_POS, largoMano);
	asigna_Array3(R_FINGERS_POS, R_WRIST_POS);
	
	add3(L_WRIST_POS, largoMano);
	asigna_Array3(L_FINGERS_POS, L_WRIST_POS);

	float largoTobillo[3] = {0.0, 0.0, HEEL_LEN};
	sub3(R_ANKLE_POS, largoTobillo);
	asigna_Array3(R_HEEL_POS, R_ANKLE_POS);
	
	sub3(L_ANKLE_POS, largoTobillo);
	asigna_Array3(L_HEEL_POS, L_ANKLE_POS);

	float largoPie[3] = {0.0, 0.0, FOOT_LEN};
	add3(R_ANKLE_POS, largoPie);
	asigna_Array3(R_TOES_POS, R_ANKLE_POS);
	
	add3(L_ANKLE_POS, largoPie);
	asigna_Array3(L_TOES_POS, R_TOES_POS);
	
	
	//////////// fin variables ragdoll ///////////////

  	int i = 0;
  	
  	dMass m;

  	dsFunctions fn;   	// setup pointers to drawstuff callback functions
  	
  	fn.version = DS_VERSION;
  	fn.start   = &start;
  	fn.step    = &simLoop;
  	fn.command = &command;
  	fn.stop    = 0;
  	fn.path_to_textures = DRAWSTUFF_TEXTURE_PATH;

  	dInitODE2(0);   	
  	
  	world        = dWorldCreate();
  	space        = dHashSpaceCreate (0);
  	contactgroup = dJointGroupCreate (0);
  	
  	dWorldSetGravity (world, 0, 0, -0.5);
  	
  	ground = dCreatePlane (space, 0, 0, 1, 0);

	cuerpos[0] = Cuerpo(world);
	body[0] = cuerpos[0].getID(); 
  	
  	cuerpos[0].setPosicion(body[0],    0,      0, STARTZ);
  		
  	dMassSetBox (&m,              1, LENGTH,  WIDTH, HEIGHT);
  	dMassAdjust (&m,          CMASS);
  	  	
  	cuerpos[0].setMasa(body[0],      &m);
  	
  	box[0] = dCreateBox (0,  LENGTH,  WIDTH, HEIGHT);
  	
  	dGeomSetBody   (box[0], body[0]);

  	for (i = 1; i < 4; i++) {   	 		
  		cuerpos[i] = Cuerpo(world);
  		body[i] = cuerpos[i].getID(); 
		
		dQuaternion q;
		dQFromAxisAndAngle (      q, 1,      0, 0, M_PI / 2);
			
		cuerpos[i].setQuaternion(body[i], q);
		dMassSetSphere (&m,          1, RADIUS);
		dMassAdjust    (&m,      WMASS);
		dBodySetMass   (body[i],    &m);
		
		sphere[i - 1] = dCreateSphere (0, RADIUS);
		
		dGeomSetBody (sphere[i - 1], body[i]);
	 }
  	
  	float lengthEntreDos = LENGTH / 2, 
  	widthEntreDos = WIDTH / 2, 
  	startZ_MenosHeightEntreDos = STARTZ - HEIGHT / 2;
  	
  	Cuerpo cuerpo1 = Cuerpo();
  		
  	cuerpo1.setPosicion(body[1],  lengthEntreDos,              0, startZ_MenosHeightEntreDos);
  	cuerpo1.setPosicion(body[2], -lengthEntreDos,  widthEntreDos, startZ_MenosHeightEntreDos);
  	cuerpo1.setPosicion(body[3], -lengthEntreDos, -widthEntreDos, startZ_MenosHeightEntreDos);
		
	Junta junta = Junta();	
		
  	for (i = 0; i < 3; i++) {   	
    		joint[i] = dJointCreateHinge2 (world, 0);
    		
    		dJointAttach (joint[i], body[0], body[i + 1]);
    		
    		const dReal *a = dBodyGetPosition(body[i + 1]);
    		    		
    		junta.setAnclaBisagra2(joint[i], a[0], a[1], a[2]); 		
    		junta.setEjesBisagra2(joint[i], zunit, yunit);
		junta.setParamsBisagra2(joint[i], dParamSuspensionERP, 0.4);		
    		junta.setParamsBisagra2(joint[i], dParamSuspensionCFM, 0.8);
    		junta.setParamsBisagra2(joint[i],dParamLoStop,0); 
    		junta.setParamsBisagra2(joint[i],dParamHiStop,0); 
	}

	Espacio espacio = Espacio(space);
	  	
	car_space = espacio.getID();  	
	
	espacio.annade(car_space,    box[0]);
	espacio.annade(car_space,    sphere[0]);
	espacio.annade(car_space,    sphere[1]);
	espacio.annade(car_space,    sphere[2]);
	
  	dMatrix3 R;
  	dRFromAxisAndAngle (       R, 0, 1,     0, -0.15);
  	
  	Caja caja = Caja(space, 2, 2, 1, R); 	//rampa  	
  	ground_box = caja.getID();

  	dsSimulationLoop (argc, argv, 800, 600, &fn);
  	
  	return 0;
}
