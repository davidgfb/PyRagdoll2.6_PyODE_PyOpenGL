class Junta {
	public:
		Junta() {
		
		}
		
		void setAnclaBisagra2(dJointID ID, dReal x, dReal y, dReal z) { 
			dJointSetHinge2Anchor (ID,  x,  y, z);
		}
		
		void setEjesBisagra2(dJointID ID, const dReal* eje, const dReal* eje1) {
			dJointSetHinge2Axes(ID,eje,eje1); 
		}
		
		void setParamsBisagra2(dJointID ID, int parameter, dReal value) {
			dJointSetHinge2Param(ID,parameter,value);
		}
};
