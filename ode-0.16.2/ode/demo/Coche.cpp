#include "Demo_Buggy.h"

class Coche {
	dJointID IDs_Juntas[3];
	Cuerpo cuerpos[3];
	dMass masa;
	dGeomID IDs_GeomsEsferas[3],
			IDs_GeomsCajas[1];
	dSpaceID ID_EspacioCoche;
	
	dReal dimensiones[3] = {LENGTH, WIDTH, HEIGHT};
	
	public: 
		Coche() {
		
		}
				
		dReal* getDims() {
			return dimensiones;
		}
			
		Coche(dWorldID ID_Mundo, dBodyID* IDs_Cuerpos, float* zunit, float* yunit, dSpaceID ID_Espacio) {
						cuerpos[0] = Cuerpo(ID_Mundo);
			IDs_Cuerpos[0] = cuerpos[0].getID(); 
		  	
		  	cuerpos[0].setPosicion(IDs_Cuerpos[0], 0, 0, STARTZ);
		  		
		  	dMassSetBox(&masa, 1, LENGTH,  WIDTH, HEIGHT);
		  	dMassAdjust(&masa, CMASS);
		  	  	
		  	cuerpos[0].setMasa(IDs_Cuerpos[0], &masa);
		  	
		  	IDs_GeomsCajas[0] = dCreateBox(0,  LENGTH,  WIDTH, HEIGHT);
		  	
		  	dGeomSetBody(IDs_GeomsCajas[0], IDs_Cuerpos[0]);

		  	for (int i = 1; i < 4; i++) {   
		  		dGeomID ID_GeomEsfera1 = dCreateSphere (0, RADIUS);
				 		
		  		cuerpos[i] = Cuerpo(ID_Mundo);
		  		IDs_Cuerpos[i] = cuerpos[i].getID(); 
				
				Cuerpo cuerpo = cuerpos[i];
				dBodyID ID_Cuerpo = IDs_Cuerpos[i];
				
				dQuaternion q;
				dQFromAxisAndAngle(q, 1, 0, 0, M_PI / 2);
					
				cuerpo.setQuaternion(IDs_Cuerpos[i], q);
				dMassSetSphere(&masa, 1, RADIUS);
				dMassAdjust(&masa, WMASS);
				dBodySetMass(ID_Cuerpo, &masa);
				
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
			
			espacio.annade(ID_EspacioCoche, IDs_GeomsCajas[0]);
			espacio.annade(ID_EspacioCoche, IDs_GeomsEsferas[0]);
			espacio.annade(ID_EspacioCoche, IDs_GeomsEsferas[1]);
			espacio.annade(ID_EspacioCoche, IDs_GeomsEsferas[2]);
		}
		
		void simula() {
		
		}
		
		dJointID getID_JuntaRuedaDelantera() {
			return IDs_Juntas[0];
		}
};
