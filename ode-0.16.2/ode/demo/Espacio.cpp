class Espacio {
	dSpaceID ID;

	public:
		Espacio(dSpaceID pID_Espacio) {
			ID = dSimpleSpaceCreate (pID_Espacio);   	// create car space and add it to the top level space
		}
		
		dSpaceID getID() {
			return ID;
		}
		
		void annade(dSpaceID ID, dGeomID ID_Geom) {
			dSpaceAdd (ID, ID_Geom);
		}
};
