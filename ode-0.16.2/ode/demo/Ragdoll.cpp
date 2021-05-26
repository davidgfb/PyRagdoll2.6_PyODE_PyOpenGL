#include "Ragdoll.h"
#include "Matriz.cpp"

class Ragdoll {
	//Creates a ragdoll of standard size at the given offset.	
	dWorldID ID_Mundo; 
	dSpaceID ID_Espacio;
	float densidad = 0.0,
	      masaTotal = 0.0,
	      offset[3] = {0.0, 0.0, 0.0}; //new float[3]
	
	public:	
		Ragdoll(dWorldID pID_Mundo, dSpaceID pID_Espacio, float pDensidad, float pOffset[3]) { //no densidad		
			ID_Mundo = pID_Mundo;
			ID_Espacio = pID_Espacio;
			densidad = pDensidad;
			asigna_Array3(offset, pOffset);			
			
			//////////// asignacion variables ////////////
	float largoAntebrazo[3] = {UPPER_ARM_LEN, 0.0, 0.0},
		largoBrazo[3] = {FORE_ARM_LEN, 0.0, 0.0},
		largoMano[3] = {HAND_LEN, 0.0, 0.0},
		largoTobillo[3] = {0.0, 0.0, HEEL_LEN},
		largoPie[3] = {0.0, 0.0, FOOT_LEN};
			
	asigna_Array3(R_ELBOW_POS, sub3(R_SHOULDER_POS, largoAntebrazo));	
	asigna_Array3(L_ELBOW_POS, add3(L_SHOULDER_POS, largoAntebrazo));	
	asigna_Array3(R_WRIST_POS, sub3(R_ELBOW_POS, largoBrazo));	
	asigna_Array3(L_WRIST_POS, add3(L_ELBOW_POS, largoBrazo));	
	asigna_Array3(R_FINGERS_POS, sub3(R_WRIST_POS, largoMano));	
	asigna_Array3(L_FINGERS_POS, add3(L_WRIST_POS, largoMano));
	asigna_Array3(R_HEEL_POS, sub3(R_ANKLE_POS, largoTobillo));	
	asigna_Array3(L_HEEL_POS, sub3(L_ANKLE_POS, largoTobillo));
	asigna_Array3(R_TOES_POS, add3(R_ANKLE_POS, largoPie));
	
	add3(L_ANKLE_POS, largoPie); //mal
	asigna_Array3(L_TOES_POS, R_TOES_POS);	
	//////////// fin asignacion variables ///////////////
		}
		
		//sobrecarga para tener parametro offset definido por defecto
		Ragdoll(dWorldID pID_Mundo, dSpaceID pID_Espacio, float pDensidad) {
			float offset[3] = {0.0, 0.0, 0.0};
			Ragdoll(ID_Mundo, pID_Espacio, pDensidad, offset);
		}
		
		void annadeCuerpo(float p_P[3], float p_P1[3], float pRadio) {
			//Adds a capsule body between joint positions p1 and p2 and with given radius to the ragdoll.
			add3(p_P, offset);
			add3(p_P1, offset);
			
			//cylinder length not including endcaps, make capsules overlap by half radius at joints
			float longCilindro = dist3(p_P, p_P1) - pRadio; //largo cilindro
		
			dBodyID ID_Cuerpo = dBodyCreate(ID_Mundo);
		
			dMass masa; 
			
			dMassSetCylinder(&masa, densidad, 3, pRadio, longCilindro); //dCreateCapsule? //3=eje z //m.
		}
};

