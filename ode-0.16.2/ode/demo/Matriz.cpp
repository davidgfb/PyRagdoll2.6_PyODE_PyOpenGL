float R_ELBOW_POS[3],
	L_ELBOW_POS[3],
	R_WRIST_POS[3],
	L_WRIST_POS[3],
	R_FINGERS_POS[3],
	L_FINGERS_POS[3],

	R_HEEL_POS[3],
	L_HEEL_POS[3],
	R_TOES_POS[3],
	L_TOES_POS[3];

bool esPositivo(float x) {
	bool positivo = false;

	if (x > 0.0) {
		positivo = true;
	}

	return positivo;
}

float sign(float x) { //sentido bool
	//Returns 1.0 if x is positive, -1.0 if x is negative or zero.
	float valor = -1.0; // (x <= 0.0), x E (-inf, 0]
	
	if (esPositivo(x)) { // (0, inf)
		valor = 1.0;
	} 
	
	return valor;
}

float len3(float v[3]) { 
	//Returns the length of 3-vector v.
	float x = v[0], y = v[1], z = v[2];
	
	return sqrt(x * x + y * y + z * z);
}

/*
void imprimeArray(float* array, int t) { //float star = tipo array	
	for (int posicion = 0; posicion < t; posicion++) { 
		printf("%f ", array[posicion]); //array);
	}	
}
*/

void asigna_Array(float *array, float *array1, int t) {
	//array=array1
	for (int pos = 0; pos < t; pos++) {
		array[pos] = array1[pos];
	}
}

float matriz3[3], 
	matriz9[9],
	matriz16[16];
		
float* neg3(float v[3]) { 
	//the negation of 3-vector v.
	float matrizNegada1[3] = {-v[0], -v[1], -v[2]};
	
	asigna_Array(matriz3, matrizNegada1,3);
			
	return matriz3;
}

float* add3(float a[3], float b[3]) { 
	//the sum of 3-vectors a and b.
	float matrizSuma1[3] = {a[0] + b[0],
				a[1] + b[1],
				a[2] + b[2]};
				
	asigna_Array(matriz3, matrizSuma1,3);
				
	return matriz3;
}

float* sub3(float a[3], float b[3]) { 
	//the difference between 3-vectors a and b.
	return add3(a, neg3(b));
}

float* mul3(float v[3], float s) { 
	//3-vector v multiplied by scalar s.
	float matrizMultiplicada1[3] = {v[0] * s, 
					v[1] * s,
					v[2] * s};
					
	asigna_Array(matriz3, matrizMultiplicada1, 3);
	
	return matriz3;
}

float* div3(float v[3], float s) {
	//3-vector v divided by scalar s.	
	return mul3(v, 1.0 / s);
}

float dist3(float a[3], float b[3]) {
	//Returns the distance between point 3-vectors a and b.		
	return len3(sub3(a, b));
}

float* norm3(float v[3]) {
	//the unit length 3-vector parallel to 3-vector v.	
	float l = len3(v);

	if (esPositivo(l)) { 
		float *matrizDividida = div3(v, l);
		
		asigna_Array(matriz3, matrizDividida, 3);
		
	} else {
		float matrizCeros[3] = {0.0, 0.0, 0.0};
		
		asigna_Array(matriz3, matrizCeros,3);
	}
	
	return matriz3;
}

float dot3(float a[3], float b[3]) {
	//Returns the dot product of 3-vectors a and b.
	float xa = a[0], ya = a[1], za = a[2], 
	      xb = b[0], yb = b[1], zb = b[2];
	      
	return xa * xb + ya * yb + za * zb;
}

float* cross(float a[3], float b[3]) { 
	//the cross product of 3-vectors a and b.
	float matrizVectorizada1[3] = {a[1] * b[2] - b[1] * a[2],
				       a[2] * b[0] - b[2] * a[0],
				       a[0] * b[1] - b[0] * a[1]};
	
	asigna_Array(matriz3, matrizVectorizada1, 3);
	
	return matriz3;
}

//float matrizProyectada[3];

float* project3(float v[3], float d[3]) {
	//projection of 3-vector v onto unit 3-vector d.
	return mul3(v, dot3(norm3(v), d));
}

float acosdot3(float a[3], float b[3]) {
	//Returns the angle between unit 3-vectors a and b.
	float x = dot3(a, b), 
		valor = 0.0; // valor E (0, inf)
	
	if (x < -1.0) { //llama
		valor = M_PI;
	} else {
		if (x == 0.0) {
			valor = acos(0.0); 
		}
	}
	
	return valor;
}

float* invert3x3(float m[9]) {
	//the inversion (transpose) of 3x3 rotation matrix m.
	float matrizInvertida1[9] = {m[0],m[3],m[6],
				      m[1],m[4],m[7],
				      m[2],m[5],m[8]};
	
	asigna_Array(matriz9, matrizInvertida1, 9);
	
	return matriz9;
}

float* rotate3(float m[9], float v[3]) { 
	//the rotation of 3-vector v by 3x3 (row major) matrix m.	
	float matrizRotada31[3] = {v[0] * m[0] + v[1] * m[1] + v[2] * m[2],
	     v[0] * m[3] + v[1] * m[4] + v[2] * m[5],
 	     v[0] * m[6] + v[1] * m[7] + v[2] * m[8]};
 	     
 	asigna_Array(matriz3, matrizRotada31, 3);     
 	     
 	return matriz3;     
}

float* zaxis(float m[9]) {
	//the z-axis vector from 3x3 (row major) rotation matrix m.
	float ejeZ1[9]; 
	
	asigna_Array(matriz9, ejeZ1, 9); //inicializa ejeZ importante!
		
	matriz9[0] = m[2];
	matriz9[1] = m[5];
	matriz9[2] = m[8];
	
	return matriz9;	
}

float* calcRotMatrix(float axis[9], float angle) {
	//the row-major 3x3 rotation matrix defining a rotation around axis by	angle. dRFromAxisAndAngle?
	float cosTheta = cos(angle), 
	      sinTheta = sin(angle),
	      t = 1.0 - cosTheta,
	      tPorEje0 = t * axis[0],
	      tPorEje1 = t * axis[1],
	      tPorEje0PorEje1 = tPorEje0 * axis[1],
	      tPorEje0PorEje2 = tPorEje0 * axis[2],
	      tPorEje1PorEje2 = tPorEje1 * axis[2],
	      senoThetaPorEje0 = sinTheta * axis[0],
	      senoThetaPorEje1 = sinTheta * axis[1],
	      senoThetaPorEje2 = sinTheta * axis[2],
	      matrizRotada91[9] = {tPorEje0 * axis[0] + cosTheta,
	      			  tPorEje0PorEje1 - senoThetaPorEje2,
	      			  tPorEje0PorEje2 + senoThetaPorEje1,
	      			  
	      			  tPorEje0PorEje1 + senoThetaPorEje2,
	      			  tPorEje1 * axis[1] + cosTheta,
	      			  tPorEje1PorEje2 - senoThetaPorEje0,
	      			  
	      			  tPorEje0PorEje2 - senoThetaPorEje1,
	      			  tPorEje1PorEje2 + senoThetaPorEje0,
	      			  t * axis[2] * axis[2] + cosTheta};
	      
	asigna_Array(matriz9,matrizRotada91,9);
	
	return matriz9;
}

float* makeOpenGLMatrix(float r[16], float p[3]) { //no puede ser float r[9] -> r[16]
	//an OpenGL compatible (column-major, 4x4 homogeneous) transformation	matrix from ODE compatible (row-major, 3x3) rotation matrix r and position	vector p.
	float s[16] = {r[0],r[3],r[6],
			0.0,r[1],r[4],
			r[7],0.0,r[2],
			r[5],r[8],0.0,
			p[0],p[1],p[2],1.0};
			
	asigna_Array(matriz16, s, 16); 
		
	return matriz16;
}

//obsoleto?
void asigna_Array3(float array[3], float array1[3]) {
	asigna_Array(array, array1, 3);
}

///////////////////// ragdoll /////////////////////
