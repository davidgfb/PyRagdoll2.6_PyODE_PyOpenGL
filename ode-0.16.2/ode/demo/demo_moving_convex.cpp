#include <ode/ode.h>
#include <drawstuff/drawstuff.h>

dWorldID ID_Mundo;
dSpaceID ID_Espacio;
dGeomID ID_GeomCapsula; //puntero a struct vacio dxGeom
dJointGroupID ID_GrupoJuntas;

////////////////// Clases //////////////////////
//////////////////// Capsula //////////////////
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
		
		/////////////// setter ///////////////		
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
//////////////////// fin Capsula //////////////////

Capsula capsula; //no se ha inicializado pero ha usado constructor vacio tiene que declararse debajo de la def clase

//////////////// Camara /////////////////////
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
			
			dsSetViewpoint(p_Pos, pRot);
			
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
//////////////// fin Camara /////////////////////

class Contacto {
	char cToString[200]; 
	dContact *arrayContactos;

	public:
		Contacto() {

		}
				
		Contacto(dContact pArrayContactos[1]) {			
			dContact arrayContactos0 = pArrayContactos[0];
			
			arrayContactos0.surface.mode = 20; //dContactBounce + dContactSoftCFM
			arrayContactos0.surface.mu = dInfinity;
			arrayContactos0.surface.mu2 = 0;
			arrayContactos0.surface.bounce = 0.1;
			arrayContactos0.surface.bounce_vel = 0.1;
			arrayContactos0.surface.soft_cfm = 0.01;
						
			arrayContactos = pArrayContactos;
			
			sprintf(cToString, "Contacto: {modoSuperficie = %i, coefFriccCoulomb = %f, coefFriccCoulombDir = %f, rebote = %f, vRebote = %f, suavidadContactoNormal = %f}\n\n", arrayContactos0.surface.mode, 
									   arrayContactos0.surface.mu, 
									   arrayContactos0.surface.mu2,
									   arrayContactos0.surface.bounce,
									   arrayContactos0.surface.bounce_vel,
									   arrayContactos0.surface.soft_cfm);
		}
				
		char* toString() {
			return cToString;
		}
};

class Cuerpo {
	char cToString[15];
	dBodyID ID;
	
	public:
		Cuerpo() {
		
		}
	
		Cuerpo(dWorldID ID_Mundo) {
			ID = dBodyCreate(ID_Mundo); 
			
			sprintf(cToString, "Cuerpo: {}\n\n");
		}
			
		dBodyID getID() {
			return ID;
		}
		
		char* toString() {
			return cToString;
		}
};
////////////////// fin Clases //////////////////////

Contacto contacto;

static void nearCallback(void *data, dGeomID o, dGeomID o1) {	
	dContact arrayContactos[1]; //es un struct  		
		
	contacto = Contacto(arrayContactos); //se reutiliza 
	//printf("%s", contacto.toString());	
		
	dCollide(o, o1, 1, &arrayContactos[0].geom, sizeof(dContact)); //requiere puntero array //&contact[0] != &contacto
			
	dJointAttach(dJointCreateContact(ID_Mundo, ID_GrupoJuntas, arrayContactos), dGeomGetBody(o), dGeomGetBody(o1));		
}

Cuerpo cuerpo;

static void start() {
	float pos[3] = {-5.0,   0.0, 5.0},
		  rot[3] = { 0.0, -45.0, 0.0}, //tilt,pan,roll
		  ladoX = 1.0, 
		  ladoY = 2.0; 	
		  	
	Camara cam = Camara(pos, rot);
	printf("%s", cam.toString());		 
		  			
	dMass m;				
	dMatrix3 R; //float a[3]
			
	cuerpo = Cuerpo(ID_Mundo);
	printf("%s", cuerpo.toString());	

	dBodyID ID_Cuerpo = cuerpo.getID();  

	dBodySetPosition(ID_Cuerpo, 0, 0, 3);	
	dRFromAxisAndAngle(R, 0.0, 0.0, 0.0, 0.0);				
	dBodySetRotation(ID_Cuerpo, R); 
	dBodySetData(ID_Cuerpo, (void*) (dsizeint) 0); 

	dMassSetCapsule(&m, 5.0, 3, ladoX, ladoY);
	
	ID_GeomCapsula = dCreateCapsule(ID_Espacio, ladoX, ladoY); 
	   
	dGeomSetBody(ID_GeomCapsula, ID_Cuerpo); 
	dBodySetMass(ID_Cuerpo, &m); 
	
	capsula = Capsula(ID_GeomCapsula); 		
	printf("%s\n", capsula.toString());
}

void drawGeom(dGeomID ID_GeomCapsula, const dReal *pos, const dReal *R, int show_aabb) {	
		dsDrawCapsule(dGeomGetPosition(ID_GeomCapsula), 					  dGeomGetRotation(ID_GeomCapsula), 					  capsula.getLargo(), 					  						  capsula.getRadio());		
}

static void simLoop(int pause) {
	dSpaceCollide(ID_Espacio, 0, &nearCallback); //sin esto no hay colisiones
	dWorldQuickStep(ID_Mundo, 0.05); //entra bucle	
	dJointGroupEmpty(ID_GrupoJuntas); //sin esto rebota y elimina g
	
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
	ID_Mundo = dWorldCreate();

	ID_Espacio = dSimpleSpaceCreate(0);
	ID_GrupoJuntas = dJointGroupCreate(0);
	dWorldSetGravity(ID_Mundo, 0, 0, -10 );
	//dWorldSetCFM(world, 1e-5);
	dCreatePlane(ID_Espacio, 0, 0, 1, 0);

	dsSimulationLoop(argc, argv, 1000, 1000, &fn); //fundamental
	
	return 0;
}

 	  	 

