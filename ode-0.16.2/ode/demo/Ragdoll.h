//rotation directions are named by the third (z-axis) row of the 3x3 matrix, because ODE capsules are oriented along the z-axis
float rightRot[9] = {0.0, 0.0, -1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0},
      leftRot[9]  = {0.0, 0.0, 1.0, 0.0, 1.0, 0.0, -1.0, 0.0, 0.0},
      upRot[9]    = {1.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 1.0, 0.0},
      downRot[9]  = {1.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 1.0, 0.0},
      bkwdRot[9]  = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0},
	//axes used to determine constrained joint rotations
	rightAxis[3] = { 1.0,  0.0,  0.0},
	leftAxis[3]  = {-1.0,  0.0,  0.0},
	upAxis[3]    = { 0.0,  1.0,  0.0},
	downAxis[3]  = { 0.0, -1.0,  0.0},
	bkwdAxis[3]  = { 0.0,  0.0,  1.0},
	fwdAxis[3]   = { 0.0,  0.0, -1.0},
	
	UPPER_ARM_LEN =  0.3,
	FORE_ARM_LEN  = 0.25,
	HAND_LEN      = 0.13, // wrist to mid-fingers only
	FOOT_LEN      = 0.18, //ankles to base of ball of foot only
	HEEL_LEN      = 0.05,
	
	BROW_H     = 1.68,
	MOUTH_H    = 1.53,
	NECK_H     =  1.5,
	SHOULDER_H = 1.37,
	CHEST_H    = 1.35,
	HIP_H      = 0.86,
	KNEE_H     = 0.48,
	ANKLE_H    = 0.08,
	
	SHOULDER_W = 0.41,
	CHEST_W    = 0.36, // actually wider, but we want narrower than shoulders (esp. with large radius)
	LEG_W      = 0.28, // between middles of upper legs
	PELVIS_W   = 0.25, // actually wider, but we want smaller than hip width;
	R_SHOULDER_POS[3] = {(float) (SHOULDER_W / -2.0), SHOULDER_H, 0.0},
	L_SHOULDER_POS[3] = {(float) (SHOULDER_W / 2.0), SHOULDER_H, 0.0},
	R_HIP_POS[3]   = {(float) (LEG_W / -2.0),   HIP_H, 0.0},
	L_HIP_POS[3]   = {(float) (LEG_W /  2.0),   HIP_H, 0.0},
	R_KNEE_POS[3]  = {(float) (LEG_W / -2.0),  KNEE_H, 0.0},
	L_KNEE_POS[3]  = {(float) (LEG_W /  2.0),  KNEE_H, 0.0},
	R_ANKLE_POS[3] = {(float) (LEG_W / -2.0), ANKLE_H, 0.0},
	L_ANKLE_POS[3] = {(float) (LEG_W /  2.0), ANKLE_H, 0.0};





























