class Coche {
	dJointID IDs_Juntas[3];

	public: 
		Coche() {
		
		}
	
		Coche(dWorldID ID_Mundo, dBodyID* IDs_Cuerpos, Junta junta, float* zunit, float* yunit) {
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
		}
		
		dJointID getID_JuntaRuedaDelantera() {
			return IDs_Juntas[0];
		}
};
