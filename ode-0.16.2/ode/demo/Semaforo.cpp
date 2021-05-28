int nFotograma = 0,
	nFotogramas = 600;
	
float posEsferaRoja[3] = {0,0,1},
	rotEsferaRoja[12] = {1.0,0.0,0.0,0.0,
			0.0,1.0,0.0,0.0,
 			0.0,0.0,1.0,0.0};

void cambiaSemaforo() {
 	if (2.0 * nFotograma < nFotogramas) { 
 		dsSetColor(1,0,0); //rojo
 	} else { //>= 		
 		dsSetColor(0,1,0); //verde	
 	}
 	
 	if (nFotograma > nFotogramas) {
 		nFotograma = 0;
 	} else {
 		 nFotograma++;
 	}
			 
  	dsDrawSphere(posEsferaRoja, rotEsferaRoja, 0.1); 
  	dsSetColor(1,1,1); //blanco (gris)
}
