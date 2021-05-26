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

static dWorldID ID_Mundo; // dynamics and collision objects (chassis, 3 wheels, environment)
static dSpaceID ID_Espacio, ID_EspacioCoche;
static dBodyID IDs_Cuerpos[4];
static dJointID IDs_Juntas[3]; 
static dJointGroupID contactgroup;
static dGeomID ground, ID_GeomCaja[1], IDs_GeomsEsferas[3], ID_GeomRampa;
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

bool esPositivo(float x) {
	bool positivo = false;

	if (x > 0.0) {
		positivo = true;
	}

	return positivo;
}


///////////// ragdoll //////////////////
float sign(float x) { //sentido bool
	//Returns 1.0 if x is positive, -1.0 if x is negative or zero.
	float valor = -1.0; // (x <= 0.0), x E (-inf, 0]
	
	if (esPositivo(x)) { // (0, inf)
		valor = 1.0;
	} 
	
	return valor;
}

float len3(float v[3]) { 
	//Returns the length of 3-vector v.
	float x = v[0], y = v[1], z = v[2];
	
	return sqrt(x * x + y * y + z * z);
}

/*
void imprimeArray(float* array, int t) { //float star = tipo array	
	for (int posicion = 0; posicion < t; posicion++) { 
		printf("%f ", array[posicion]); //array);
	}	
}
*/

void asigna_Array(float *array, float *array1, int t) {
	//array=array1
	for (int pos = 0; pos < t; pos++) {
		array[pos] = array1[pos];
	}
}

float matriz3[3], 
	matriz9[9],
	matriz16[16];
		
float* neg3(float v[3]) { //no devuelve nada haz una copia antes
	//the negation of 3-vector v.
	float matrizNegada1[3] = {-v[0], -v[1], -v[2]};
	
	asigna_Array(matriz3,
			matrizNegada1,3);
			
	return matriz3;
}

float* add3(float a[3], float b[3]) { //modifica a haz una copia
	//the sum of 3-vectors a and b.
	float matrizSuma1[3] = {a[0] + b[0],
				a[1] + b[1],
				a[2] + b[2]};
				
	asigna_Array(matriz3,
			matrizSuma1,3);
				
	return matriz3;
}

float* sub3(float a[3], float b[3]) { //modifica a haz una copia
	//the difference between 3-vectors a and b.
	return add3(a, neg3(b));
}

float* mul3(float v[3], float s) { 
	//3-vector v multiplied by scalar s.
	float matrizMultiplicada1[3] = {v[0] * s, 
					v[1] * s,
					v[2] * s};
					
	asigna_Array(matriz3,
			matrizMultiplicada1, 3);
	
	return matriz3;
}

float* div3(float v[3], float s) {
	//3-vector v divided by scalar s.	
	return mul3(v, 1.0 / s);
}

float dist3(float a[3], float b[3]) {
	//Returns the distance between point 3-vectors a and b.		
	return len3(sub3(a, b));
}

float* norm3(float v[3]) {
	//the unit length 3-vector parallel to 3-vector v.	
	float l = len3(v);

	if (esPositivo(l)) { 
		float *matrizDividida = div3(v, l);
		
		asigna_Array(matriz3,
			matrizDividida,3);
		
	} else {
		float matrizCeros[3] = {0.0, 0.0, 0.0};
		
		asigna_Array(matriz3,
			matrizCeros,3);
	}
	
	return matriz3;
}

float dot3(float a[3], float b[3]) {
	//Returns the dot product of 3-vectors a and b.
	float xa = a[0], ya = a[1], za = a[2], 
	      xb = b[0], yb = b[1], zb = b[2];
	      
	return xa * xb + ya * yb + za * zb;
}

float* cross(float a[3], float b[3]) { 
	//the cross product of 3-vectors a and b.
	float matrizVectorizada1[3] = {a[1] * b[2] - b[1] * a[2],
				       a[2] * b[0] - b[2] * a[0],
				       a[0] * b[1] - b[0] * a[1]};
	
	asigna_Array(matriz3,
		     matrizVectorizada1,3);
	
	return matriz3;
}

//float matrizProyectada[3];

float* project3(float v[3], float d[3]) {
	//projection of 3-vector v onto unit 3-vector d.
	return mul3(v, dot3(norm3(v), d));
}

float acosdot3(float a[3], float b[3]) {
	//Returns the angle between unit 3-vectors a and b.
	float x = dot3(a, b), valor = 0.0; // valor E (0, inf)
	
	if (x < -1.0) { //llama
		valor = M_PI;
	} else {
		if (x == 0.0) {
			valor = acos(0.0); 
		}
	}
	
	return valor;
}

float* invert3x3(float m[9]) {
	//the inversion (transpose) of 3x3 rotation matrix m.
	float matrizInvertida1[9] = {m[0],m[3],m[6],
				      m[1],m[4],m[7],
				      m[2],m[5],m[8]};
	
	asigna_Array(matriz9,
		     matrizInvertida1,9);
	
	return matriz9;
}

float* rotate3(float m[9], float v[3]) { 
	//the rotation of 3-vector v by 3x3 (row major) matrix m.	
	float matrizRotada31[3] = {v[0] * m[0] + v[1] * m[1] + v[2] * m[2],
	     v[0] * m[3] + v[1] * m[4] + v[2] * m[5],
 	     v[0] * m[6] + v[1] * m[7] + v[2] * m[8]};
 	     
 	asigna_Array(matriz3,matrizRotada31,3);     
 	     
 	return matriz3;     
}

float* zaxis(float m[9]) {
	//the z-axis vector from 3x3 (row major) rotation matrix m.
	float ejeZ1[9]; 
	
	asigna_Array(matriz9,ejeZ1,9); //inicializa ejeZ importante!
		
	matriz9[0] = m[2];
	matriz9[1] = m[5];
	matriz9[2] = m[8];
	
	return matriz9;	
}

float* calcRotMatrix(float axis[9], float angle) {
	//the row-major 3x3 rotation matrix defining a rotation around axis by	angle. dRFromAxisAndAngle?
	float cosTheta = cos(angle), 
	      sinTheta = sin(angle),
	      t = 1.0 - cosTheta,
	      tPorEje0 = t * axis[0],
	      tPorEje1 = t * axis[1],
	      tPorEje0PorEje1 = tPorEje0 * axis[1],
	      tPorEje0PorEje2 = tPorEje0 * axis[2],
	      tPorEje1PorEje2 = tPorEje1 * axis[2],
	      senoThetaPorEje0 = sinTheta * axis[0],
	      senoThetaPorEje1 = sinTheta * axis[1],
	      senoThetaPorEje2 = sinTheta * axis[2],
	      matrizRotada91[9] = {tPorEje0 * axis[0] + cosTheta,
	      			  tPorEje0PorEje1 - senoThetaPorEje2,
	      			  tPorEje0PorEje2 + senoThetaPorEje1,
	      			  
	      			  tPorEje0PorEje1 + senoThetaPorEje2,
	      			  tPorEje1 * axis[1] + cosTheta,
	      			  tPorEje1PorEje2 - senoThetaPorEje0,
	      			  
	      			  tPorEje0PorEje2 - senoThetaPorEje1,
	      			  tPorEje1PorEje2 + senoThetaPorEje0,
	      			  t * axis[2] * axis[2] + cosTheta};
	      
	asigna_Array(matriz9,matrizRotada91,9);
	
	return matriz9;
}

float* makeOpenGLMatrix(float r[16], float p[3]) { //no puede ser float r[9] -> r[16]
	//an OpenGL compatible (column-major, 4x4 homogeneous) transformation	matrix from ODE compatible (row-major, 3x3) rotation matrix r and position	vector p.
	float s[16] = {r[0],r[3],r[6],
			0.0,r[1],r[4],
			r[7],0.0,r[2],
			r[5],r[8],0.0,
			p[0],p[1],p[2],1.0};
			
	asigna_Array(matriz16,s,16); 
		
	return matriz16;
}

//obsoleto?
void asigna_Array3(float array[3], float array1[3]) {
	asigna_Array(array,array1,3);
}

//////////////////////////////////////

dJointID getJuntaRuedaDelantera() {
	return IDs_Juntas[0];
}

static void nearCallback (void *, dGeomID o1, dGeomID o2) {
	  int i = 0, n = 0;

	  bool g1 = (o1 == ground || o1 == ID_GeomRampa), 	  // only collide things with the ground
	  g2 = (o2 == ground || o2 == ID_GeomRampa);
	  
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
			      
			      dJointID c = dJointCreateContact (ID_Mundo, contactgroup, &contact[i]);
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
}

static void command (int cmd) {  // lower 
					
}

int ciclosEspera = 20, cicloActual = 0;

////////////////// simLoop ////////////////////////
static void simLoop (int pause) { // simulation loop
   	int i = 0;
	float vMax = 0.1;
		
	dJointID juntaRuedaDelanteraCoche = getJuntaRuedaDelantera();	
		  
    	if (!pause) {
         	dJointSetHinge2Param (juntaRuedaDelanteraCoche,  dParamVel2, -speed); // motor rueda/s delantera
         	dJointSetHinge2Param (juntaRuedaDelanteraCoche, dParamFMax2,    0.1);

    	 	dReal v = steer - dJointGetHinge2Angle1 (getJuntaRuedaDelantera()); // steering
     	 	
    	 	dJointSetHinge2Param (getJuntaRuedaDelantera(),         dParamVel,     v);
    	 	dJointSetHinge2Param (getJuntaRuedaDelantera(),        dParamFMax,   0.2);
    	 	dJointSetHinge2Param (getJuntaRuedaDelantera(),      dParamLoStop, -0.75);
    	 	dJointSetHinge2Param (getJuntaRuedaDelantera(),      dParamHiStop,  0.75);
    	 	dJointSetHinge2Param (juntaRuedaDelanteraCoche, dParamFudgeFactor,   0.1);

    	 	dSpaceCollide(ID_Espacio,    0, &nearCallback);
    	 	dWorldStep(ID_Mundo, 0.05);

    	 	dJointGroupEmpty(contactgroup); // remove all contact joints
    	}

  	dReal sides[3] = {LENGTH, WIDTH, HEIGHT};
  	
  	dsDrawBox (dBodyGetPosition(IDs_Cuerpos[0]), dBodyGetRotation(IDs_Cuerpos[0]), sides);
  		
	for (i = 1; i < 4; i++) {
  		 dsDrawCylinder(dBodyGetPosition(IDs_Cuerpos[i]),
				 dBodyGetRotation(IDs_Cuerpos[i]), 0.02, RADIUS);
	}

  	dVector3 ss;
  	dGeomBoxGetLengths (ID_GeomRampa, ss);
  	
  	dsDrawBox (dGeomGetPosition(ID_GeomRampa), dGeomGetRotation(ID_GeomRampa), ss);
}
////////////////// fin simLoop ////////////////////

///////////////////// ragdoll /////////////////////


class Ragdoll {
	//Creates a ragdoll of standard size at the given offset.
	//constructor
	//public: 	
	
	dWorldID ID_Mundo; 
	dSpaceID ID_Espacio;
	float densidad = 0.0,

	      totalMass = 0.0,
		  offset[3] = {0.0, 0.0, 0.0}; //new float[3]
	
	
	Ragdoll(dWorldID w, dSpaceID s, float d, float o[3]) { //no densidad		
		ID_Mundo = w;
		ID_Espacio = s;
		densidad = d;
		asigna_Array3(offset, o); //offset = {o[0],o[1],o[2]};
		
	}
	
	//sobrecarga para tener parametro offset definido por defecto
	Ragdoll(dWorldID ID_Mundo, dSpaceID pID_Espacio, float density) {
		float offset[3] = {0.0, 0.0, 0.0};
		Ragdoll(ID_Mundo,pID_Espacio,density,offset);
	}
	
	void addBody(float p1[3], float p2[3], float radius) {
		//Adds a capsule body between joint positions p1 and p2 and with given radius to the ragdoll.
		add3(p1, offset);
		add3(p2, offset);
		
		//cylinder length not including endcaps, make capsules overlap by half radius at joints
		float cyllen = dist3(p1, p2) - radius; //largo cilindro
	
		dBodyID body = dBodyCreate(ID_Mundo);
	
		dMass m; 
		
		dMassSetCylinder(&m,densidad,3,radius,cyllen); //dCreateCapsule? //3=eje z //m.
	}
};


class Caja {
	dGeomID ID;

	public:
		Caja(dSpaceID pID_Espacio, dReal dx, dReal dy, dReal dz, dMatrix3 R) {
			ID = dCreateBox(pID_Espacio, dx, dy, dz);   
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
		Espacio(dSpaceID pID_Espacio) {
			ID = dSimpleSpaceCreate (pID_Espacio);   	// create car space and add it to the top level space
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

  	dsFunctions llamadas_Simulacion;
  	
  	llamadas_Simulacion.version = DS_VERSION;
  	llamadas_Simulacion.start   = &start;
  	llamadas_Simulacion.step    = &simLoop;
  	llamadas_Simulacion.command = &command;
  	llamadas_Simulacion.stop    = 0;
  	llamadas_Simulacion.path_to_textures = DRAWSTUFF_TEXTURE_PATH;

  	dInitODE2(0);   	
  	
  	ID_Mundo = dWorldCreate();
  	ID_Espacio = dHashSpaceCreate(0);
  	contactgroup = dJointGroupCreate(0);
  	
  	dWorldSetGravity(ID_Mundo, 0, 0, -0.5);
  	
  	ground = dCreatePlane(ID_Espacio, 0, 0, 1, 0);

	cuerpos[0] = Cuerpo(ID_Mundo);
	IDs_Cuerpos[0] = cuerpos[0].getID(); 
  	
  	cuerpos[0].setPosicion(IDs_Cuerpos[0], 0, 0, STARTZ);
  		
  	dMassSetBox(&m, 1, LENGTH,  WIDTH, HEIGHT);
  	dMassAdjust(&m, CMASS);
  	  	
  	cuerpos[0].setMasa(IDs_Cuerpos[0], &m);
  	
  	ID_GeomCaja[0] = dCreateBox(0,  LENGTH,  WIDTH, HEIGHT);
  	
  	dGeomSetBody(ID_GeomCaja[0], IDs_Cuerpos[0]);

  	for (i = 1; i < 4; i++) {   
  		dGeomID ID_GeomEsfera1 = dCreateSphere (0, RADIUS);
		 		
  		cuerpos[i] = Cuerpo(ID_Mundo);
  		IDs_Cuerpos[i] = cuerpos[i].getID(); 
		
		Cuerpo cuerpo = cuerpos[i];
		dBodyID ID_Cuerpo = IDs_Cuerpos[i];
		
		dQuaternion q;
		dQFromAxisAndAngle(q, 1, 0, 0, M_PI / 2);
			
		cuerpo.setQuaternion(IDs_Cuerpos[i], q);
		dMassSetSphere(&m, 1, RADIUS);
		dMassAdjust(&m, WMASS);
		dBodySetMass(ID_Cuerpo, &m);
		
		cuerpos[i] = cuerpo;
		
		IDs_GeomsEsferas[i - 1] = ID_GeomEsfera1;
		
		dGeomSetBody(ID_GeomEsfera1, ID_Cuerpo);
	 }
  	
  	float lengthEntreDos = LENGTH / 2, 
	  	widthEntreDos = WIDTH / 2, 
	  	startZ_MenosHeightEntreDos = STARTZ - HEIGHT / 2;
  	
  	Cuerpo cuerpo = Cuerpo();
  		
  	cuerpo.setPosicion(IDs_Cuerpos[1],  lengthEntreDos,              0, startZ_MenosHeightEntreDos);
  	cuerpo.setPosicion(IDs_Cuerpos[2], -lengthEntreDos,  widthEntreDos, startZ_MenosHeightEntreDos);
  	cuerpo.setPosicion(IDs_Cuerpos[3], -lengthEntreDos, -widthEntreDos, startZ_MenosHeightEntreDos);
		
	Junta junta = Junta();	
		
  	for (i = 0; i < 3; i++) {   
  		dJointID ID_Junta = dJointCreateHinge2 (ID_Mundo, 0);
  		
    		IDs_Juntas[i] = ID_Junta;
    		
    		dJointAttach (ID_Junta, IDs_Cuerpos[0], IDs_Cuerpos[i + 1]);
    		
    		const dReal *a = dBodyGetPosition(IDs_Cuerpos[i + 1]);
    		    		
    		junta.setAnclaBisagra2(ID_Junta, a[0], a[1], a[2]); 		
    		junta.setEjesBisagra2(ID_Junta, zunit, yunit);
		junta.setParamsBisagra2(ID_Junta, dParamSuspensionERP, 0.4);		
    		junta.setParamsBisagra2(ID_Junta, dParamSuspensionCFM, 0.8);
    		junta.setParamsBisagra2(ID_Junta,dParamLoStop,0); 
    		junta.setParamsBisagra2(ID_Junta,dParamHiStop,0); 
	}

	Espacio espacio = Espacio(ID_Espacio);
	  	
	ID_EspacioCoche = espacio.getID();  	
	
	espacio.annade(ID_EspacioCoche, ID_GeomCaja[0]);
	espacio.annade(ID_EspacioCoche, IDs_GeomsEsferas[0]);
	espacio.annade(ID_EspacioCoche, IDs_GeomsEsferas[1]);
	espacio.annade(ID_EspacioCoche, IDs_GeomsEsferas[2]);
	
  	dMatrix3 R;
  	dRFromAxisAndAngle (R, 0, 1, 0, -0.15);
  	
  	Caja caja = Caja(ID_Espacio, 2, 2, 1, R); 	//rampa  	
  	ID_GeomRampa = caja.getID();

  	dsSimulationLoop (argc, argv, 800, 600, &llamadas_Simulacion);
  	
  	return 0;
}
