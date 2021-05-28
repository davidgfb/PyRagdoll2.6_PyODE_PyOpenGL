/*
buggy with suspension.
this also shows you how to use geom groups.
*/

#include <ode/ode.h>
#include <drawstuff/drawstuff.h>
#include "Ragdoll.cpp"
#include "Caja.cpp"
#include "Espacio.cpp"
#include "Cuerpo.cpp"
#include "Junta.cpp"
#include "Coche.cpp"
#include "Semaforo.cpp"

dVector3 yunit = {0, 1, 0}, 
			zunit = {0, 0, 1};
dWorldID ID_Mundo; 

dBodyID IDs_Cuerpos[4];
dJointGroupID ID_GrupoJunta;
dGeomID ID_GeomSuelo, 	 
	ID_GeomRampa;
dSpaceID ID_Espacio;
float dxCaja = 2,
	  	dyCaja = dxCaja,
	  	dzCaja = 1,
	  	
	  	xCaja = dxCaja,
	  	yCaja = 0,
	  	zCaja = -0.34,
	  	
	  	posCaja[3] = {xCaja, yCaja, zCaja},
	  	dimsCaja[3] = {dxCaja, dyCaja, dzCaja};
dMatrix3 rotCaja;

Coche coche;

static void nearCallback (void *, dGeomID o1, dGeomID o2) {
	  bool g1 = (o1 == ID_GeomSuelo || o1 == ID_GeomRampa), 	  // only collide things with the ground
	  		g2 = (o2 == ID_GeomSuelo || o2 == ID_GeomRampa);
	  
	  if (g1 ^ g2) {
		  dContact contactos[10];
		  int nContactos = dCollide (o1, o2, 10, &contactos[0].geom, sizeof(dContact));
		  
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
		
float speed = 1.5;
			
void simLoop (int pause) { 	
	dJointID juntaRuedaDelanteraCoche = coche.getID_JuntaRuedaDelantera();	
		
	dJointSetHinge2Param (juntaRuedaDelanteraCoche,dParamVel2,-speed);
    dJointSetHinge2Param (juntaRuedaDelanteraCoche,dParamFMax2,0.1);
	
		/*  
         	dJointSetHinge2Param(juntaRuedaDelanteraCoche, dParamFMax2, 0.1);    	 	
    	 	dJointSetHinge2Param(juntaRuedaDelanteraCoche, dParamFMax, 0.2);
    	 	dJointSetHinge2Param(juntaRuedaDelanteraCoche, dParamLoStop, -0.75);
    	 	dJointSetHinge2Param(juntaRuedaDelanteraCoche, dParamHiStop, 0.75);
    	 	dJointSetHinge2Param(juntaRuedaDelanteraCoche, dParamFudgeFactor, 0.1);
		*/
		
    	 	dSpaceCollide(ID_Espacio, 0, &nearCallback);
    	 	dWorldStep(ID_Mundo, 0.05); //10 g?

    	 	dJointGroupEmpty(ID_GrupoJunta); //


  	dReal sides[3] = {LENGTH, WIDTH, HEIGHT}; 	
  	dBodyID ID_Cuerpo0 = IDs_Cuerpos[0];
   	  
   	cambiaSemaforo();
  	
  	dsDrawBox(dBodyGetPosition(ID_Cuerpo0), dBodyGetRotation(ID_Cuerpo0), sides);
  	
	for (int i = 0; i < 4; i++) {
		dBodyID ID_Cuerpo = IDs_Cuerpos[i];
  		dsDrawCylinder(dBodyGetPosition(ID_Cuerpo),
			       dBodyGetRotation(ID_Cuerpo), 
			       0.02, RADIUS);
	}
  	
  	dsDrawBox(posCaja, rotCaja, dimsCaja);
}

//clase Moto...

int main(int argc, char **argv) {  	
  	dsFunctions llamadas_Simulacion;
  	
  	llamadas_Simulacion.version = DS_VERSION;
  	llamadas_Simulacion.start   = &start;
  	llamadas_Simulacion.step    = &simLoop;
  	llamadas_Simulacion.command = 0;
  	llamadas_Simulacion.stop    = 0;
  	llamadas_Simulacion.path_to_textures = 0; 

  	dInitODE2(0);   	
  	
  	ID_Mundo = dWorldCreate();
  	ID_Espacio = dHashSpaceCreate(0);
  	ID_GrupoJunta = dJointGroupCreate(0);
  	
  	dWorldSetGravity(ID_Mundo, 0, 0, -0.5); //
  	
  	ID_GeomSuelo = dCreatePlane(ID_Espacio, 0, 0, 1, 0);
			
  	coche = Coche(ID_Mundo,IDs_Cuerpos,zunit,yunit,ID_Espacio);

  	dRFromAxisAndAngle (rotCaja, 0, 1, 0, -0.15);
  		
  	Caja rampa = Caja(ID_Espacio, dxCaja, dyCaja, dzCaja, xCaja, yCaja, zCaja, rotCaja); 	  	
  	ID_GeomRampa = rampa.getID();

  	dsSimulationLoop(argc, argv, 800, 600, &llamadas_Simulacion);
  	  	
  	return 0;
}
