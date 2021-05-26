class Caja {
	dGeomID ID;

	public:
		Caja(dSpaceID pID_Espacio, dReal dx, dReal dy, dReal dz, dMatrix3 R) {
			ID = dCreateBox(pID_Espacio, dx, dy, dz);   
			dGeomSetPosition(ID, 2, 0, -0.34);
  			dGeomSetRotation(ID, R);
		}
		
		dGeomID getID() {
			return ID;
		}
};
