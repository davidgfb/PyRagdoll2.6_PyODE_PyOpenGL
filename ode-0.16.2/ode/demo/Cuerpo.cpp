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
