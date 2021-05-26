/*
buggy with suspension.
this also shows you how to use geom groups.
*/
#include <ode/ode.h>
#include <drawstuff/drawstuff.h>
#include "texturepath.h"
#include "Ragdoll.h"
#include "Ragdoll.cpp"
#include "Caja.cpp"
#include "Espacio.cpp"
#include "Cuerpo.cpp"
#include "Junta.cpp"

float   LENGTH = 0.7, 
	WIDTH = 0.5, 
	HEIGHT = 0.2, 
	RADIUS = 0.18, 
	STARTZ = 0.5, 
	CMASS = 1.0, 
	WMASS = 0.2; // some constants

dVector3 yunit = {0, 1, 0}, 
			zunit = {0, 0, 1};

dWorldID ID_Mundo; // dynamics and collision objects (chassis, 3 wheels, environment)
dSpaceID ID_Espacio, ID_EspacioCoche;
dBodyID IDs_Cuerpos[4];
dJointID IDs_Juntas[3]; 
dJointGroupID ID_GrupoJunta;
dGeomID ground, ID_GeomsCajas[1], IDs_GeomsEsferas[3], ID_GeomRampa;
dReal speed = 0, 
		steer = 0; // things that the user controls // user commands

Cuerpo cuerpos[3];

dJointID getJuntaRuedaDelantera() {
	return IDs_Juntas[0];
}

static void nearCallback (void *, dGeomID o1, dGeomID o2) {
	  int nContactos = 0;
	  bool g1 = (o1 == ground || o1 == ID_GeomRampa), 	  // only collide things with the ground
	  		g2 = (o2 == ground || o2 == ID_GeomRampa);
	  
	  if (g1 ^ g2) {
		  dContact contactos[10];
		  nContactos = dCollide (o1, o2, 10, &contactos[0].geom, sizeof(dContact));
		  
		  int sumaContactos = dContactSlip1 + dContactSlip2 +
				      dContactSoftERP + dContactSoftCFM + dContactApprox1;
		  
			    for (int i = 0; nContactos > 0 && i < nContactos; i++) {
			    	dContact contacto = contactos[i];
			    
			      contacto.surface.mode = sumaContactos;
			      
			      contacto.surface.mu = dInfinity;
			      
			      contacto.surface.slip1    = 0.1;
			      contacto.surface.slip2    = 0.1;
			      contacto.surface.soft_erp = 0.5;
			      contacto.surface.soft_cfm = 0.3;
			      
			      dJointID c = dJointCreateContact (ID_Mundo, ID_GrupoJunta, &contactos[i]); //
			      dJointAttach (c, dGeomGetBody(contacto.geom.g1),
					       dGeomGetBody(contacto.geom.g2));
			    }

	  }
}

void start() { 
	  float pos[3] = {1, -1, 1}, 
	               rot[3] = {121.0, -27.5, 0.0};
	  
	  dsSetViewpoint(pos, rot);
}

void simLoop (int pause) { // simulation loop		
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

    	 	dJointGroupEmpty(ID_GrupoJunta); //
    	}

  	dReal sides[3] = {LENGTH, WIDTH, HEIGHT};
  	
  	dBodyID ID_Cuerpo0 = IDs_Cuerpos[0];
  	
  	dsDrawBox(dBodyGetPosition(ID_Cuerpo0), dBodyGetRotation(ID_Cuerpo0), sides);
  		
	for (int i = 1; i < 4; i++) {
		dBodyID ID_Cuerpo = IDs_Cuerpos[i];
  		 dsDrawCylinder(dBodyGetPosition(ID_Cuerpo),
				 dBodyGetRotation(ID_Cuerpo), 0.02, RADIUS);
	}

  	dVector3 dimensionesLateralesCaja;
  	dGeomBoxGetLengths(ID_GeomRampa, dimensionesLateralesCaja);
  	
  	dsDrawBox (dGeomGetPosition(ID_GeomRampa), dGeomGetRotation(ID_GeomRampa), dimensionesLateralesCaja);
}

//clase Coche, Moto...

int main(int argc, char **argv) {  	
  	dMass m;

  	dsFunctions llamadas_Simulacion;
  	
  	llamadas_Simulacion.version = DS_VERSION;
  	llamadas_Simulacion.start   = &start;
  	llamadas_Simulacion.step    = &simLoop;
  	llamadas_Simulacion.command = 0;
  	llamadas_Simulacion.stop    = 0;
  	llamadas_Simulacion.path_to_textures = DRAWSTUFF_TEXTURE_PATH;

  	dInitODE2(0);   	
  	
  	ID_Mundo = dWorldCreate();
  	ID_Espacio = dHashSpaceCreate(0);
  	ID_GrupoJunta = dJointGroupCreate(0);
  	
  	dWorldSetGravity(ID_Mundo, 0, 0, -0.5); //
  	
  	ground = dCreatePlane(ID_Espacio, 0, 0, 1, 0);

	cuerpos[0] = Cuerpo(ID_Mundo);
	IDs_Cuerpos[0] = cuerpos[0].getID(); 
  	
  	cuerpos[0].setPosicion(IDs_Cuerpos[0], 0, 0, STARTZ);
  		
  	dMassSetBox(&m, 1, LENGTH,  WIDTH, HEIGHT);
  	dMassAdjust(&m, CMASS);
  	  	
  	cuerpos[0].setMasa(IDs_Cuerpos[0], &m);
  	
  	ID_GeomsCajas[0] = dCreateBox(0,  LENGTH,  WIDTH, HEIGHT);
  	
  	dGeomSetBody(ID_GeomsCajas[0], IDs_Cuerpos[0]);

  	for (int i = 1; i < 4; i++) {   
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
		
  	for (int i = 0; i < 3; i++) {   
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
	
	espacio.annade(ID_EspacioCoche, ID_GeomsCajas[0]);
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
