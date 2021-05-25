#include <ode/ode.h>
#include <drawstuff/drawstuff.h>

dWorldID ID_Mundo;
dSpaceID ID_Espacio;
dGeomID ID_GeomCapsula; //puntero a struct vacio dxGeom
dJointGroupID ID_GrupoJuntas;

////////////////// Clases //////////////////////
class Capsula { 
	char cToString[100]; 
	float radioCapsula,
		  largoCapsula;
	dGeomID ID_Geometrico;
	
	void setParamsGeoms(dGeomID pID_GeomCapsula) {									
			dGeomCapsuleGetParams(pID_GeomCapsula, &radioCapsula, &largoCapsula); //Setter?
	}
	
	public:	
		//////////// constructor ////////////////////
		Capsula() { //para declarar clase no inicializada
			
		}
	
		Capsula(dSpaceID pID_Espacio, float pRadio, float pLargo) {
			ID_Geometrico = dCreateCapsule(pID_Espacio, pRadio, pLargo); 
			setParamsGeoms(ID_Geometrico);			
			
			sprintf(cToString, "Capsula: {radio = %f, largo = %f}\n\n", getRadio(), getLargo());
		}
		////////////// fin constructor ////////////
		
		void dibuja() {
			dsDrawCapsule(dGeomGetPosition(ID_GeomCapsula), 				  	      dGeomGetRotation(ID_GeomCapsula), 				          getLargo(), 					  					  	      getRadio());
		}
		
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
		
		dGeomID getID_Geometrico() {
			return ID_Geometrico;
		}
		//////////// fin getter ///////////////
			
		char* toString() {			
			return cToString;
		}
};

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
				
		Contacto(dContact pArrayContactos[1], int pModoSup, float pCoefFriccCoulomb, float pCoefFriccCoulombDir, float pRebote, float pV_Rebote, float pSuavidadContactoNormal) {			
			dContact arrayContactos0 = pArrayContactos[0];
			
			arrayContactos0.surface.mode = pModoSup; //dContactBounce + dContactSoftCFM
			arrayContactos0.surface.mu = pCoefFriccCoulomb;
			arrayContactos0.surface.mu2 = pCoefFriccCoulomb;
			arrayContactos0.surface.bounce = pRebote;
			arrayContactos0.surface.bounce_vel = pV_Rebote;
			arrayContactos0.surface.soft_cfm = pSuavidadContactoNormal; //suavidadContactoNormal = pSuavidadContactoNormal
						
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
	
		Cuerpo(dWorldID pID_Mundo, dMatrix3 R, dReal x, dReal y, dReal z) {
			ID = dBodyCreate(pID_Mundo); 
			
			dBodySetPosition(ID, x, y, z);		
			dBodySetRotation(ID, R); 
			dBodySetData(ID, (void*) 0); //puntero vacio a 0
			
			sprintf(cToString, "Cuerpo: {}\n\n");
		}
			
		dBodyID getID() {
			return ID;
		}
		
		char* toString() {
			return cToString;
		}
};

class Objeto {

};

class Matriz {
	char cToString[20];
	dMatrix3 R;

	public:		
		Matriz() {				
			sprintf(cToString, "Matriz: {}\n\n");
		}
		
		void rotaRadianesEnEje(dReal x, dReal y, dReal z, dReal anguloRadianes) {
			dRFromAxisAndAngle(R, x, y, z, anguloRadianes);
		}
		
		dReal* getMatriz() {
			return R;
		}
		
		char* toString() {
			return cToString;
		}
};

class Masa {
	char cToString[20];
	dMass masa;

	public:
		Masa() {
		
		}
	
		Masa(dBodyID pID_Cuerpo, dReal pDensidad, int pDireccion, float pLadoX, float pLadoY) {
			dMassSetCapsule(&masa, pDensidad, pDireccion, pLadoX, pLadoY);
			dBodySetMass(pID_Cuerpo, &masa); 
						
			sprintf(cToString, "Masa: {}\n\n");
		}
				
		char* toString() {
			return cToString;
		}
};

class Geometria {
	char cToString[20];

	public:
		Geometria(dGeomID pID_GeomCapsula, dBodyID pID_Cuerpo, dSpaceID ID_Espacio, float ladoX, float ladoY) {					
			dGeomSetBody(pID_GeomCapsula, pID_Cuerpo); 
		
			sprintf(cToString, "Geometria: {}\n\n");
		}
				
		char* toString() {
			return cToString;
		}
};
////////////////// fin Clases //////////////////////

Contacto contacto;

static void nearCallback(void *data, dGeomID o, dGeomID o1) {	
	dContact arrayContactos[1]; //es un struct  		
		
	contacto = Contacto(arrayContactos, 20, dInfinity, 0.0, 0.1, 0.1, 0.01); //se reutiliza 
	//printf("%s", contacto.toString());	
		
	dCollide(o, o1, 1, &arrayContactos[0].geom, sizeof(dContact)); //requiere puntero array //&contact[0] != &contacto
			
	dJointAttach(dJointCreateContact(ID_Mundo, ID_GrupoJuntas, arrayContactos), dGeomGetBody(o), dGeomGetBody(o1));		
}

Cuerpo cuerpo;

static void start() {
	float pos[3] = {-5.0,   0.0, 5.0},
		  rot[3] = { 0.0, -45.0, 0.0}, //tilt,pan,roll
		  radioCapsula = 1.0, 
		  largoCapsula = 2.0; 	
		  
	dReal xCuerpo = 0, yCuerpo = 0, zCuerpo = 3;
		  	
	Camara cam = Camara(pos, rot);
			  		  					
	Matriz matriz = Matriz();	  
		
	matriz.rotaRadianesEnEje(0.0, 0.0, 0.0, 0.0);
				
	cuerpo = Cuerpo(ID_Mundo, matriz.getMatriz(), xCuerpo, 
												  yCuerpo, 
												  zCuerpo);

	dBodyID ID_Cuerpo = cuerpo.getID();  			
					   
	capsula = Capsula(ID_Espacio, radioCapsula, 
								  largoCapsula); 
	/*
	radioCapsula = capsula.getRadio(),
	largoCapsula = capsula.getLargo();  
	*/
								  
	Masa masa = Masa(ID_Cuerpo, 5.0, 3, radioCapsula, 
										largoCapsula);
	   
	ID_GeomCapsula = capsula.getID_Geometrico();	   
	   
	Geometria geometria = Geometria(ID_GeomCapsula, ID_Cuerpo, ID_Espacio, radioCapsula, largoCapsula);
	   	
	printf("%s%s%s%s%s%s", cam.toString(), 
						   matriz.toString(), 
						   cuerpo.toString(),
						   masa.toString(), 
						   geometria.toString(), 
						   capsula.toString());
}

static void simLoop(int pause) {
	dSpaceCollide(ID_Espacio, 0, &nearCallback); //sin esto no hay colisiones
	dWorldQuickStep(ID_Mundo, 0.05); //entra bucle	
	dJointGroupEmpty(ID_GrupoJuntas); //sin esto rebota y elimina g
	
	capsula.dibuja();
}

class Plano {
	public:
		Plano(dSpaceID espacio, dReal a, dReal b, dReal c, dReal d) {
			dCreatePlane(ID_Espacio, 0, 0, 1, 0);
		}
};

class Mundo {
	dWorldID ID;
	
	public:
		Mundo() {
			ID = dWorldCreate();
		}
		
		dWorldID getID() {
			return ID;
		}
};

class Espacio {
	dSpaceID ID;
	
	public:
		Espacio() {
			ID = dSimpleSpaceCreate(0);
		}
		
		dSpaceID getID() {
			return ID;
		}
};

int main(int argc, char **argv) {
	dsFunctions fn;
	
	fn.version = DS_VERSION;
	fn.start = &start;
	fn.step = &simLoop;
	fn.stop = 0;
	fn.path_to_textures = DRAWSTUFF_TEXTURE_PATH;

	dInitODE2(0);
	
	Mundo mundo = Mundo();
	
	ID_Mundo = mundo.getID();

	Espacio espacio = Espacio();
	ID_Espacio = espacio.getID(); //objeto.getID() padre
		
	ID_GrupoJuntas = dJointGroupCreate(0);
	dWorldSetGravity(ID_Mundo, 0, 0, -10 );
	//dWorldSetCFM(world, 1e-5);
	
	Plano plano = Plano(ID_Espacio, 0, 0, 1, 0);

	dsSimulationLoop(argc, argv, 1000, 1000, &fn); //fundamental
	
	return 0;
}

 	  	 

