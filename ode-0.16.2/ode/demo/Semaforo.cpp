class Semaforo {
	int nFotograma = 0,
		nFotogramas = 600;
	
	float posEsfera[3] = {0,0, 2},
		rotEsfera[12] = {1.0, 0.0,0.0,0.0,
			0.0, 1.0, 0.0,0.0,
 			0.0,0.0, 1.0, 0.0};
 	public:		
	 	Semaforo() {
	 	
	 	}

		//void cambiaColor() {
		
		//}
		
		void actualiza() { //cuentaFotogramas
		 	if (2.0 * nFotograma < nFotogramas) { 
		 		dsSetColor(1, 0,0); //rojo
		 	} else { //>= 		
		 		dsSetColor(0, 1, 0); //verde	
		 	}
		 	
		 	if (nFotograma > nFotogramas) {
		 		nFotograma = 0;
		 	} else { //<=
		 		 nFotograma++;
		 	}
					 
		  	dsDrawSphere(posEsfera, rotEsfera, 0.1); 
		  	dsSetColor(1,1,1); //blanco (gris)
		}
};
