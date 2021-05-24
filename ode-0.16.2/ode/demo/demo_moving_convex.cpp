#include <ode/ode.h>
#include <drawstuff/drawstuff.h>

dWorldID world;
dSpaceID space;
dGeomID ID_GeomCapsula; //puntero a struct vacio dxGeom
dJointGroupID contactgroup;

static void nearCallback(void *data, dGeomID o, dGeomID o1) {
	dContact contact[1], //es un struct  		
		     contacto = contact[0];
		
	contacto.surface.mode = 20; //dContactBounce + dContactSoftCFM
	contacto.surface.mu = dInfinity;
	contacto.surface.mu2 = 0;
	contacto.surface.bounce = 0.1;
	contacto.surface.bounce_vel = 0.1;
	contacto.surface.soft_cfm = 0.01;
		
	dCollide(o, o1, 1, &contact[0].geom, sizeof(dContact));
		//&contact[0] != &contacto
	dRSetIdentity(new dMatrix3);
			
	dJointAttach(dJointCreateContact(world, contactgroup, contact), dGeomGetBody(o), dGeomGetBody(o1));		
}

class Capsula { 
	char cToString[100]; 
	float radioCapsula,
		  largoCapsula;
	
	void setParamsGeoms(dGeomID pID_GeomCapsula) {									
			dGeomCapsuleGetParams(pID_GeomCapsula, &radioCapsula, &largoCapsula); //Setter?
	}
	
	public:	
		//////////// constructor ////////////////////
		Capsula() { //para declarar clase no inicializada
			
		}
	
		Capsula(dGeomID pID_GeomCapsula) {
			setParamsGeoms(pID_GeomCapsula);
			
			sprintf(cToString, "Capsula: {radio = %f, largo = %f}\n", getRadio(), getLargo());
		}
		////////////// fin constructor ////////////
		
		/////////////// setter /////////////////		
		void setID_GeomCapsula(dGeomID pID_GeomCapsula) {
			ID_GeomCapsula = pID_GeomCapsula;
		}
		///////////// fin setter //////////////
		
		////////////// getter /////////////////
		float getRadio() {
			return radioCapsula;			
		}
		
		float getLargo() {
			return largoCapsula;
		}
		//////////// fin getter ///////////////
			
		char* toString() {			
			return cToString;
		}
};

Capsula capsula; //no se ha inicializado pero ha usado constructor vacio tiene que declararse debajo de la def clase

class Camara {
	float pos[3], //unidades x,y,z
		  rot[3]; //angulos alfa,beta,gamma
	char cToString[200]; 
	
	public:
		Camara(float p_Pos[3], float pRot[3]) {
			for (int posicion = 0; posicion < 3; posicion++) {
				pos[posicion] = p_Pos[posicion];
				rot[posicion] = pRot[posicion];
			} 	
			
			sprintf(cToString, "Camara: {x = %f, y = %f, z = %f,\nalfa = %f, beta = %f, gamma = %f}\n\n", 
					getX(), getY(), getZ(),
					getAlfa(), getBeta(), getGamma());	
		}
		
		/////////////// getter ///////////////////		
		///////////// posicion /////////////////
		float getX() {
			return pos[0];
		}
		
		float getY() {
			return pos[1];
		}
		
		float getZ() {
			return pos[2];
		}
		///////////// fin posicion ///////////////		
		///////////// rotacion //////////////////
		float getAlfa() {
			return rot[0];
		}
		
		float getBeta() {
			return rot[1];
		}
		
		float getGamma() {
			return rot[2];
		}
		///////////// fin rotacion ////////////////
		//////////////// fin getter /////////////
		
		char* toString() {			
			return cToString;
		}
};

/*
char* arrayFloats_A_ArrayChars(float arrayFloats[]) {
	char arrayChars[50];
  	int n = sprintf(arrayChars, "%f", arrayFloats); //len(arrayChars)
  	
	return arrayChars;
}
*/

static void start() {
	float xyz[3] = {-5.0,   0.0, 5.0},
		  tpr[3] = { 0.0, -45.0, 0.0}, //tilt,pan,roll
		  ladoX = 1.0, 
		  ladoY = 2.0; 	
		  
	Camara cam = Camara(xyz, tpr);
	printf("%s", cam.toString());		
	
	/*
	printf("%sx = %f, y = %f, z = %f}\n", cam.toString(), cam.getX(), cam.getY(), cam.getZ());
	*/
	/*
	for (int pos = 0; pos < 3; pos++) {
		printf("%f", );
	}
	*/
	//printf("%s\n", (char*) cam.getPos());		  
		  			
	dMass m;				
	dMatrix3 R; //float a[3]
	dBodyID obj = dBodyCreate(world);  	

	dsSetViewpoint(xyz, tpr);
	dBodySetPosition(obj, 0, 0, 3);	
	dRFromAxisAndAngle(R, 0.0, 0.0, 0.0, 0.0);				
	dBodySetRotation(obj, R); 
	dBodySetData(obj, (void*) (dsizeint) 0); 

	dMassSetCapsule(&m, 5.0, 3, ladoX, ladoY);
	
	ID_GeomCapsula = dCreateCapsule(space, ladoX, ladoY); 
	   
	dGeomSetBody(ID_GeomCapsula, obj); 
	dBodySetMass(obj, &m); 
	
	capsula = Capsula(ID_GeomCapsula); //ID_GeomCapsula
		
	printf("%s\n", capsula.toString());
	//printf("%p\n", capsula.getID_GeomCapsula());
}

void drawGeom(dGeomID ID_GeomCapsula, const dReal *pos, const dReal *R, int show_aabb) {	
		dsDrawCapsule(dGeomGetPosition(ID_GeomCapsula), 					  dGeomGetRotation(ID_GeomCapsula), 					  capsula.getLargo(), 					  						  capsula.getRadio());		
}

static void simLoop(int pause) {
	dSpaceCollide(space, 0, &nearCallback); //sin esto no hay colisiones
	dWorldQuickStep(world, 0.05); //entra bucle	
	dJointGroupEmpty(contactgroup); //sin esto rebota y elimina g
	
	drawGeom(ID_GeomCapsula, 0, 0, 0); //pinta 
}

int main(int argc, char **argv) {
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
	dWorldSetGravity(world, 0, 0, -10 );
	//dWorldSetCFM(world, 1e-5);
	dCreatePlane(space, 0, 0, 1, 0);

	dsSimulationLoop(argc, argv, 1000, 1000, &fn); //fundamental
	
	return 0;
}

 	  	 

