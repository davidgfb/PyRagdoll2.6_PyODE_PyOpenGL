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
	dMatrix3 RI;
	dRSetIdentity(RI);
			
	dJointAttach(dJointCreateContact(world, contactgroup, contact), dGeomGetBody(o), dGeomGetBody(o1));		
}

class Capsula { 
	float parametros[2]; //private
	//char txt[] = "{Capsula: "+(char*) ID_GeomCapsula;
	char sCapsula[10] = "{Capsula:";	
	dGeomID ID_GeomCapsula;
	//char aC_ID_GeomCapsula[] = (char*) ID_GeomCapsula;
	
	public:	
		Capsula() { 
			//para declarar clase no inicializada
		}
	
		Capsula(dGeomID pID_GeomCapsula) {
			setID_GeomCapsula(pID_GeomCapsula);
		}
		
		////////////// getter /////////////////
		dGeomID getID_GeomCapsula() {
			return ID_GeomCapsula;
		}
		//////////// fin getter ///////////////
		
		/////////////// setter /////////////////
		void setID_GeomCapsula(dGeomID pID_GeomCapsula) {
			ID_GeomCapsula = pID_GeomCapsula;
		}
		///////////// fin setter //////////////
		
		float* getParamsGeoms() {
			//devuelve puntero a array float2 de parametros geoms capsula
			float radio_Obtenido, 
				  largo_Obtenido;
			
			dGeomCapsuleGetParams(getID_GeomCapsula(), &radio_Obtenido, &largo_Obtenido);

			parametros[0] = radio_Obtenido;
			parametros[1] = largo_Obtenido;

			return parametros; 
		}
		
		float getRadio() {
			float *parametros = getParamsGeoms(),
				  radioCapsula = parametros[0];
				  
			return radioCapsula;			
		}
		
		float getLargo() {
			float *parametros = getParamsGeoms(),
				  largoCapsula = parametros[1];
		
			return largoCapsula;
		}
			
		char* toString() {			
			return sCapsula;
		}
};

Capsula capsula; //no se ha inicializado pero ha usado constructor vacio

static void start() {
	float xyz[3] = {-5.0,   0.0, 5.0},
		  tpr[3] = { 0.0, -45.0, 0.0},
		  ladoX = 1.0, 
		  ladoY = 2.0; //tilt,pan,roll				
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

 	  	 

