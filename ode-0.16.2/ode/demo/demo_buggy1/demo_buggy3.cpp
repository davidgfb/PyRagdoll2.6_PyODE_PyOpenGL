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

#ifdef dDOUBLE // select correct drawing functions
	#define dsDrawBox dsDrawBoxD
	#define dsDrawSphere dsDrawSphereD
	#define dsDrawCylinder dsDrawCylinderD
	#define dsDrawCapsule dsDrawCapsuleD
#endif

//////////// variables ////////////////////
float LENGTH = 0.7, WIDTH = 0.5, HEIGHT = 0.2, RADIUS = 0.18, STARTZ = 0.5, CMASS = 1.0, WMASS = 0.2; // some constants

static const dVector3 yunit = {0, 1, 0}, zunit = {0, 0, 1};

static dWorldID world; // dynamics and collision objects (chassis, 3 wheels, environment)
static dSpaceID space, car_space;
static dBodyID body[4];

static dJointID joint[3]; // array3d de dJointIDs

static dJointGroupID contactgroup;

static dGeomID ground, box[1], sphere[3], ground_box;

static dReal speed = 0, steer = 0; // things that the user controls // user commands
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
bool esPositivo(float x) {
	//Devuelve verdadero si x es positivo, falso si es negativo o cero.
	bool positivo = false; // (x <= 0.0), x E (-inf, 0]
	
	if (x > 0.0) { // (0, inf)
		positivo = true;
	}
	
	return positivo;
}

//obsoleto
float sign(float x) { 
	//Devuelve 1.0 si x es positivo, -1.0 si es negativo o cero.
	float signo = -1.0; 
	
	if (esPositivo(x)) { 
		signo = 1.0;
	} 
	
	return signo;
}

float sentido(float x) { //signo 
	return sign(x);
}

float len3(float v[3]) { 
	//Returns the length of 3-vector v.
	float x = v[0], y = v[1], z = v[2];
	
	return sqrt(x * x + y * y + z * z);
}

//////////////// procedimiento /////////////////
//error sizeof(puntero) // tama??o
void imprimeArray(float* array, int t) { //float star = tipo array	
	for (int posicion = 0; posicion < t; posicion++) { 
		printf("%f ", array[posicion]); //array);
	}	
}
/////////////// fin procedimiento /////////////////

void asigna_ArrayT(float* a, float* b, int t) {
	for (int posicion = 0; posicion < t; posicion++) {
		a[posicion] = b[posicion];
	}
	/*
	array[0] = array1[0];
	array[1] = array1[1];
	array[2] = array1[2];
	*/
}

//obsoleta
void asigna_Array3(float a[3], float b[3]) {
	asigna_ArrayT(a,b,3);
}

//obsoleta
void asigna_Array(float a[3], float b[3]) {
	asigna_Array3(a,b);
}

float v1[]; //cuidado

float* neg(float* v, int t) {
	//the negation of vector v. 
	asigna_ArrayT(v1,v,t); //hace una copia
	
	for (int posicion = 0; posicion < 3; posicion++) {
		v1[posicion] *= -1.0; 
	}
	
	return v1;
}

//float v1[3] = {0.0, 0.0, 0.0};

float* neg3(float v[3]) { 
	return neg(v,3);
	
	//the negation of 3-vector v.
	/*
	//v1 = {v[0],v[1],v[2]}; //hace una copia
	asigna_ArrayT(&v1[0], &v[0], 3);
	
	for (int posicion = 0; posicion < 3; posicion++) {
		v1[posicion] *= -1.0; 
	}
	
	return v1;
	*/
}

float* add3(float a[3], float b[3]) { //modifica a haz una copia
	//the sum of 3-vectors a and b.
	a[0] += b[0];
	a[1] += b[1];
	a[2] += b[2];
	
	return a;
}

float* sub3(float a[3], float b[3]) { //modifica a haz una copia
	//the difference between 3-vectors a and b.
	float* bNegada = neg3(b);
	float* resultado = add3(a, bNegada);
	
	/*
	a[0] -= b[0];
	a[1] -= b[1];
	a[2] -= b[2];
	*/
	return resultado;
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



//////////////////////////////////////

dJointID getJuntaRuedaDelantera() {
	return joint[0];
}

static void nearCallback (void *, dGeomID o1, dGeomID o2) {
	  /*
	  // this is called by dSpaceCollide when two objects in space are
	  // potentially colliding.
	  */
	  int i = 0, n = 0;

	  bool g1 = (o1 == ground || o1 == ground_box); 	  // only collide things with the ground
	  bool g2 = (o2 == ground || o2 == ground_box);
	  
	  if (g1 ^ g2) {
		  const int N = 10;
		  dContact contact[N];
		  n = dCollide (o1, o2, N, &contact[0].geom, sizeof(dContact));
		  
		  int sumaContactos = dContactSlip1   + dContactSlip2   +
				      dContactSoftERP + dContactSoftCFM + dContactApprox1;
		  
		  if (n > 0) {
			    for (i = 0; i < n; i++) {
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

bool estaConduciendoCoche = false;

float vAbs(float v) {
	if (v < 0) {
		v *= -1.0;
	}
	
	return v;
}

float escalarV_R = 0.3; // velocidad a la que rueda > vMin

bool estaParadoCoche(float v) { //[0]
	bool estaParado = false;
	
	if (v == 0) {
		estaParado = true;
	} 
	
	return estaParado;
}

bool estaV_MinCoche(float v) { //vAnormalmenteReducida
	bool estaV_Min = false;
	
	if (not estaParadoCoche(v) and vAbs(v) < escalarV_R) {
		estaV_Min = true; //frena
	}
	
	return estaV_Min;
}

bool estaRodandoCoche(float v) {
	bool estaRodando = false;
	
	if (vAbs(v) == escalarV_R) { //[vR]
		estaRodando = true;
	}
	
	return estaRodando;
}

//parado = 0, vAnormalmenteReducida = (abs(+-(0)), abs(+-vR)) (frenara solo), rodando = abs(+-vR), en marcha = [abs(+-vR), abs(+-inf))

static void command (int cmd) {  // lower 
	estaConduciendoCoche = true;
		
	if (cmd == 'w' || cmd == 'W') { // flechas //acelerador
		if (estaParadoCoche(speed)) { // si esta parado en parking mete directa y empieza a rodar // puede estar rodando en neutral...
			speed = escalarV_R; // mete directa y rueda hacia delante [vR]
		} else {
			if (estaRodandoCoche(speed)) { // si esta en marcha atras o en directa acelera
				speed *= 1.3; // +-30%
			} 
		}
	}
			
	if (cmd == 's' || cmd == 'S') { //freno / marcha atras
		if (estaParadoCoche(speed) || speed >= escalarV_R) { //directa +
			speed -= escalarV_R; // mete y rueda marcha atras   //frena hacia atras		
		} else {
			if (speed <= -escalarV_R) { //marcha atras // -
				speed += escalarV_R; //frena hacia delante
			}
		}		
	}
	
	if (cmd == ' ') { // espacio = freno de mano
    		speed = 0; // bloquea las ruedas (el coche no tiene abs)
  	}
	
	/////////// giro ///////////////
	if (cmd == 'a' || cmd == 'A') { 
		if (steer > -1.0) {
	    		steer -= 0.5;
    		}
	}
	
	if (cmd == 'd' || cmd == 'D') {
		if (steer < 1.0) {
    			steer += 0.5;
    		}
	}
	////////////////////////////////// 	
}

int ciclosEspera = 20, cicloActual = 0;

////////////////// simLoop ////////////////////////
static void simLoop (int pause) { // simulation loop
   	int i = 0;
	float vMax = 0.1;
	
	/////////// freno automatico ///////////////
	if (estaV_MinCoche(speed)) { //velocidad anormalmente reducida (0, vR)
		//tiraFrenoMano() //echaFrenoMano()
		speed = 0;
	}
	////////////////////////////////////////////
	
	/////// funcion retardada corrige direccion ///////// 
	if (cicloActual > ciclosEspera) { // por si nos hemos saltado algun ciclo	
		if (estaConduciendoCoche) { //jugador.estaConduciendo(coche)
			estaConduciendoCoche = false;
		} else { //no esta conduciendo coche 
			steer = 0;
		}
		
		cicloActual = 0;
	} else {
		cicloActual++;
	}
	/////////////////////////////////////////////////
	  
    	if (!pause) {
         	dJointSetHinge2Param (getJuntaRuedaDelantera(),  dParamVel2, -speed); // motor rueda/s delantera
         	dJointSetHinge2Param (getJuntaRuedaDelantera(), dParamFMax2,    0.1);

    	 	dReal v = steer - dJointGetHinge2Angle1 (getJuntaRuedaDelantera()); // steering
    
    	 	if (v > vMax) {
    		 	v = vMax;
   	 	}
    
    	 	if (v < -vMax) { 
    	 	 	v = -vMax;
    	 	}
    
    	 	v *= 10.0;
    	 	
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
	//arrayBodies
	//arrayGeoms
	//arrayJuntas
	      totalMass = 0.0,
		  offset[3] = {0.0, 0.0, 0.0}; //new float[3]
	
	
	Ragdoll(dWorldID w, dSpaceID s, float d, float o[3]) { //no densidad		
		world = w;
		space = s;
		densidad = d;
		asigna_Array(offset, o); //offset = {o[0],o[1],o[2]};
		
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
	
		dMass m; 
		
		dMassSetCapsule(&m,densidad,3,radius,cyllen); //3=eje z 
		
		dBodySetMass(body, &m);
		
		// set parameters for drawing the body
		//dGeom dCreateCapsule dGeomCapsuleSetParams (dGeomID
		
		// create a capsule geom for collision detection
		//// draw a geom
		//void drawGeom( dGeomID g, const dReal *pos, const dReal *R, int show_aabb )
			//dGeomCapsuleGetParams( g,&radius,&length );
			
		// define body rotation automatically from body axis
		sub3(p2,p1);
		norm3(p2);
		float ceros[3] = {0.0,0.0,0.0}, //
		      za[3] = {0.0,0.0,0.0},
		      xa[3] = {0.0,1.0,0.0},
		      ejeX[3] = {1.0,0.0,0.0},
		      ya[3] = {0.0,0.0,0.0}; //?
		      
		asigna_Array(za,p2);
				
		float pEscalar = dot3(za, ejeX);
		
		if (abs(pEscalar) < 0.7) {
			asigna_Array(xa,ejeX); 
		} 
		
		//asigna_Array(ya,);
	}
};




/////////////// main ///////////////////////
int main (int argc, char **argv) {
	//////////// variables ragdoll ////////////
	float largoAntebrazo[3] = {UPPER_ARM_LEN, 0.0, 0.0};
	sub3(R_SHOULDER_POS, largoAntebrazo);
	asigna_Array(R_ELBOW_POS, R_SHOULDER_POS);

	add3(L_SHOULDER_POS, largoAntebrazo);	
	asigna_Array(L_ELBOW_POS, L_SHOULDER_POS);
	
	float largoBrazo[3] = {FORE_ARM_LEN, 0.0, 0.0};
	sub3(R_ELBOW_POS, largoBrazo);
	asigna_Array(R_WRIST_POS, R_ELBOW_POS);
	
	add3(L_ELBOW_POS, largoBrazo);
	asigna_Array(L_WRIST_POS, L_ELBOW_POS);
	
	float largoMano[3] = {HAND_LEN, 0.0, 0.0};
	sub3(R_WRIST_POS, largoMano);
	asigna_Array(R_FINGERS_POS, R_WRIST_POS);
	
	add3(L_WRIST_POS, largoMano);
	asigna_Array(L_FINGERS_POS, L_WRIST_POS);

	float largoTobillo[3] = {0.0, 0.0, HEEL_LEN};
	sub3(R_ANKLE_POS, largoTobillo);
	asigna_Array(R_HEEL_POS, R_ANKLE_POS);
	
	sub3(L_ANKLE_POS, largoTobillo);
	asigna_Array(L_HEEL_POS, L_ANKLE_POS);

	float largoPie[3] = {0.0, 0.0, FOOT_LEN};
	add3(R_ANKLE_POS, largoPie);
	asigna_Array(R_TOES_POS, R_ANKLE_POS);
	
	add3(L_ANKLE_POS, largoPie);
	asigna_Array(L_TOES_POS, R_TOES_POS);
	
	
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

  	dInitODE2(0);   	// create world
  	
  	world        = dWorldCreate();
  	space        = dHashSpaceCreate (0);
  	contactgroup = dJointGroupCreate (0);
  	
  	dWorldSetGravity (world, 0, 0, -0.5);
  	
  	ground = dCreatePlane (space, 0, 0, 1, 0);

  	body[0] = dBodyCreate (world);   	// chassis body
  	
  	dBodySetPosition (body[0],    0,      0, STARTZ);
  	dMassSetBox (&m,              1, LENGTH,  WIDTH, HEIGHT);
  	dMassAdjust (&m,          CMASS);
  	dBodySetMass  (body[0],      &m);
  	
  	box[0] = dCreateBox (0,  LENGTH,  WIDTH, HEIGHT);
  	
  	dGeomSetBody   (box[0], body[0]);

  	for (i = 1; i < 4; i++) {   	// wheel bodies
    		body[i] = dBodyCreate (world);
		
		dQuaternion q;
		dQFromAxisAndAngle (      q, 1,      0, 0, M_PI / 2);
		dBodySetQuaternion (body[i], q);
		dMassSetSphere (&m,          1, RADIUS);
		dMassAdjust    (&m,      WMASS);
		dBodySetMass   (body[i],    &m);
		
		sphere[i - 1] = dCreateSphere (0, RADIUS);
		
		dGeomSetBody (sphere[i - 1], body[i]);
	 }
  	
  	float lengthEntreDos = LENGTH / 2, widthEntreDos = WIDTH / 2, startZ_MenosHeightEntreDos = STARTZ - HEIGHT / 2;
  	
  	dBodySetPosition (body[1],  lengthEntreDos,              0, startZ_MenosHeightEntreDos);
  	dBodySetPosition (body[2], -lengthEntreDos,  widthEntreDos, startZ_MenosHeightEntreDos);
  	dBodySetPosition (body[3], -lengthEntreDos, -widthEntreDos, startZ_MenosHeightEntreDos);

  	for (i = 0; i < 3; i++) {   	// front and back wheel hinges
    		joint[i] = dJointCreateHinge2 (world, 0);
    		
    		dJointAttach (joint[i], body[0], body[i + 1]);
    		
    		const dReal *a = dBodyGetPosition (body[i + 1]);
    		
    		dJointSetHinge2Anchor (joint[i],  a[0],  a[1], a[2]);
    		dJointSetHinge2Axes   (joint[i], zunit, yunit);
  	}

  	for (i = 0; i < 3; i++) {   	// set joint suspension
    		dJointSetHinge2Param (joint[i], dParamSuspensionERP, 0.4);
    		dJointSetHinge2Param (joint[i], dParamSuspensionCFM, 0.8);
  	}

  	for (i = 1; i < 3; i++) {   	// lock back wheels along the steering axis
    		dJointSetHinge2Param (joint[i],dParamLoStop,0); //,dParamVel,0); // set stops to make sure wheels always stay in alignment
    		dJointSetHinge2Param (joint[i],dParamHiStop,0); // ,dParamFMax,dInfinity);
    		/*
    		// the following alternative method is no good as the wheels may get out
    		// of alignment:
    		*/
	}

  	car_space = dSimpleSpaceCreate (space);   	// create car space and add it to the top level space
  	dSpaceSetCleanup (car_space, 0);
  	
  	dSpaceAdd (car_space,    box[0]);
  	dSpaceAdd (car_space, sphere[0]);
  	dSpaceAdd (car_space, sphere[1]);
  	dSpaceAdd (car_space, sphere[2]);

  	ground_box = dCreateBox (space, 2, 1.5, 1);   	// environment
  	dMatrix3 R;
  	dRFromAxisAndAngle (       R, 0, 1,     0, -0.15);
  	dGeomSetPosition (ground_box, 2, 0, -0.34);
  	dGeomSetRotation (ground_box, R);

  	int w = 800, h = 600;   	// run simulation
  	dsSimulationLoop (argc, argv, w, h, &fn);

  	dGeomDestroy (   box[0]);
  	dGeomDestroy (sphere[0]);
  	dGeomDestroy (sphere[1]);
  	dGeomDestroy (sphere[2]);
  	
  	dJointGroupDestroy (contactgroup);
  	dSpaceDestroy (space);
  	dWorldDestroy (world);
  	dCloseODE();
  	
  	return 0;
}
