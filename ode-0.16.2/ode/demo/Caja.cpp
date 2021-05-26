class Caja {
	dGeomID ID;

	public:
		Caja(dSpaceID pID_Espacio, dReal dx, dReal dy, dReal dz, dReal x, dReal y, dReal z, dMatrix3 R) {
			ID = dCreateBox(pID_Espacio, dx, dy, dz);   
			dGeomSetPosition(ID, x, y, z);
  			dGeomSetRotation(ID, R);
		}
		
		dGeomID getID() {
			return ID;
		}
};
