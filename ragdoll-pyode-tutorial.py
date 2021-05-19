from time import time,sleep
from random import uniform
from sys import exit
from math import pi,sqrt,acos,cos,sin

from OpenGL.GL import glClearColor,glClear,\
GL_COLOR_BUFFER_BIT,GL_DEPTH_BUFFER_BIT,glEnable,\
GL_DEPTH_TEST,GL_LIGHTING,GL_NORMALIZE,glShadeModel,\
GL_SMOOTH,glMatrixMode,GL_PROJECTION,glLoadIdentity,\
glViewport,GL_MODELVIEW,glLightfv,GL_LIGHT0,\
GL_POSITION,GL_DIFFUSE,GL_SPECULAR,GL_COLOR_MATERIAL,\
glColor3f,glPushMatrix,glMultMatrixd,glBegin,\
GL_QUAD_STRIP,glNormal3f,glVertex3f,glEnd,\
glTranslated,glPopMatrix,\
\
glRotatef,glFlush

from OpenGL.GLU import gluPerspective,gluLookAt

from OpenGL.GLUT import glutInit,glutInitDisplayMode,\
GLUT_RGB,GLUT_DEPTH,GLUT_DOUBLE,\
glutInitWindowPosition,glutInitWindowSize,\
glutCreateWindow,glutKeyboardFunc,glutDisplayFunc,\
glutIdleFunc,glutMainLoop,glutPostRedisplay,\
glutSolidSphere,glutSwapBuffers,\
\
glutMouseFunc,glutPassiveMotionFunc

from ode import Infinity,World,Space,GeomPlane,\
JointGroup,Body,Mass,GeomCCylinder,FixedJoint,\
BallJoint,UniversalJoint,ParamLoStop,ParamHiStop,\
ParamLoStop2,ParamHiStop2,HingeJoint,areConnected,\
collide,ContactJoint

def sign(x):
	"""
	Returns 1.0 if x is positive, 
	-1.0 if x is negative or zero.
	"""  
	valor = -1.0
	
	if x > 0.0:
		valor = 1.0
	   
	return valor

def len3(v):
	"""
	Returns the length of 3-vector v.
	"""
	suma = 0.0

	for nColumna in range(3):
		suma += v[nColumna] ** 2

	return sqrt(suma)

def neg3(v):
	"""
	Returns the negation of 3-vector v.
	"""
	vector = []

	for nColumna in range(3):
		vector.append(-v[nColumna])

	return vector

def add3(a, b):
	"""
	Returns the sum of 3-vectors a and b.
	"""
	vector = []

	for nColumna in range(3):
		vector.append(a[nColumna] + b[nColumna])

	return vector

def sub3(a, b):
	"""
	Returns the difference between 3-vectors a and b.
	"""
	vector = []

	for nColumna in range(3):
		vector.append(a[nColumna] - b[nColumna])

	return vector

def mul3(v, s):
	"""
	Returns 3-vector v multiplied by scalar s.
	"""
	vector = []

	for nColumna in range(3):
		vector.append(s * v[nColumna])

	return vector

def div3(v, s):
	"""
	Returns 3-vector v divided by scalar s.
	"""
	vector = []

	for nColumna in range(3):
		vector.append(v[nColumna] / s) #no es lo mismo que s/v

	return vector

def dist3(a, b):
	"""
	Returns the distance between point 3-vectors a 
	and b.
	"""
	return len3(sub3(a, b))

def norm3(v):
	"""
	Returns the unit length 3-vector parallel to 
	3-vector v.
	"""
	l = len3(v)
	vector = [0.0, 0.0, 0.0]

	if (l > 0.0): 
		vector = div3(v,l)
	
	return vector

def dot3(a, b):
	"""
	Returns the dot product of 3-vectors a and b.
	"""
	suma = 0.0

	for nColumna in range(3):
		suma += a[nColumna] * b[nColumna]

	return suma

def cross(a, b):
	"""
	Returns the cross product of 3-vectors a and b.
	return numpy.cross(a,b)
	"""
	matriz = []

	matrizColumnas = [[1,2], #1%3,2%3
					  [2,0], #5%3,6%3  #2%3,3%3
					  [0,1]  #9%3,10%3 #3%3,1%3
					 ]

	for nFila in range(3):
		columnas = matrizColumnas[nFila]
		columna0 = columnas[0] #o matrizColumnas[nFila][0]
		columna1 = columnas[1] #o matrizColumnas[nFila][1]
		matriz.append(a[columna0] * b[columna1] - \
					  b[columna0] * a[columna1]) 

	return matriz 
	
def project3(v, d):
	"""
	Returns projection of 3-vector v onto unit 
	3-vector d.
	"""
	return mul3(v, dot3(norm3(v), d))

def acosdot3(a, b):
	"""
	Returns the angle between unit 3-vectors a and b.
	o
	valor = acos(x) # mayor carga mejor legibilidad
	elif x > 1.0:
		valor = 0.0
	"""
	x = dot3(a, b)
	valor = acos(x) # -1.0 <= x <= 1.0 [-1, 1]

	if x < -1.0: # (-inf, -1)
		valor = pi

	elif x > 1.0: # (1, inf)
		valor = 0.0

	return valor

def rotate3(m, v):
	"""
	Returns the rotation of 3-vector v by 3x3 
	(row major) matrix m.
	
	matriz = []
	matrizColumnas = [[0,0], # 0
					  [1,1],
					  [2,2],

					  [0,3], # 1,3
					  [1,4],
					  [2,5],

					  [0,6], # 2,6
					  [1,7],
					  [2,8]
					 ]

	for nGrupo in range(3): # nGrupo = nFilaMatriz 
		#nGrupo *= 3
		nSubgrupo = 3 * nGrupo
		grupo = matrizColumnas[nSubgrupo : \
							   nSubgrupo + 2]
		suma = 0.0

		for subgrupo in grupo:
			suma += v[subgrupo[0]] * m[subgrupo[1]]

		matriz.append(suma)

	return matriz
	"""
	return (v[0] * m[0] + \
			v[1] * m[1] + \
			v[2] * m[2], \
			
			v[0] * m[3] + \
			v[1] * m[4] + \
			v[2] * m[5], \
			
			v[0] * m[6] + \
			v[1] * m[7] + \
			v[2] * m[8])

def invert3x3(m):
	"""
	Returns the inversion (transpose) of 3x3 rotation 
	matrix m.
	"""
	sColumnas = '036147258'
	matriz = []

	for sColumna in sColumnas:
		posColumna = int(sColumna)
		matriz.append(m[posColumna])

	return matriz

def zaxis(m):
	"""
	Returns the z-axis vector from 3x3 (row major) 
	rotation matrix m.
	"""
	matriz = []
	sColumnas = '258'

	for sColumna in sColumnas:
		posColumna = int(sColumna)
		matriz.append(m[posColumna])

	return matriz

def calcRotMatrix(axis, angle):
	"""
	Returns the row-major 3x3 rotation matrix 
	defining a rotation around axis by
	angle.
	"""
	cosTheta = cos(angle)
	sinTheta = sin(angle)
	t = 1.0 - cosTheta
	
	return (t * axis[0] ** 2 + cosTheta,
		
		t * axis[0] * axis[1] - sinTheta * axis[2],
		t * axis[0] * axis[2] + sinTheta * axis[1],
		t * axis[0] * axis[1] + sinTheta * axis[2],
		
		t * axis[1] ** 2 + cosTheta,
		
		t * axis[1] * axis[2] - sinTheta * axis[0],
		t * axis[0] * axis[2] - sinTheta * axis[1],
		t * axis[1] * axis[2] + sinTheta * axis[0],
		
		t * axis[2] ** 2 + cosTheta)

def makeOpenGLMatrix(r, p):
	"""
	Returns an OpenGL compatible (column-major, 4x4 
	homogeneous) transformation
	matrix from ODE compatible (row-major, 3x3) 
	rotation matrix r and position
	vector p.
	
	matrizColumnas = [[0,3,6, 0.0],
					  [1,4,7, 0.0],
					  [2,5,8, 0.0],

					  [0,1,2, 1.0]
					 ]
	matriz = []

	for nFila in range(3):
		fila = matrizColumnas[nFila] 

		for nColumna in range(3):			
			matriz.append(r[fila[nColumna]])

		matriz.append(fila[3])

	for nColumna in range(3):
		matriz.append(p[3][nColumna]) 

	matriz.append(matrizColumnas[3][3]) # p[0], p[1], p[2], 1.0)
	"""

	return (r[0], r[3], r[6], 0.0,
			r[1], r[4], r[7], 0.0,
			r[2], r[5], r[8], 0.0,
			
			p[0], p[1], p[2], 1.0)
	
def getBodyRelVec(b, v):
	"""
	Returns the 3-vector v transformed into the local 
	coordinate system of ODE	body b.
	"""
	return rotate3(invert3x3(b.getRotation()), v)

'''
rotation directions are named by the third (z-axis) 
row of the 3x3 matrix,
because ODE capsules are oriented along the z-axis
HARDCODEADO!!!!!!!!
'''
rightRot        = (0.0, 0.0, -1.0, 0.0, 1.0, 0.0, 1.0,\
				   0.0, 0.0)
leftRot         = (0.0, 0.0, 1.0, 0.0, 1.0, 0.0, -1.0,\
				   0.0, 0.0)
upRot = downRot = (1.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, \
				   1.0, 0.0)
bkwdRot         = (1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, \
				   0.0, 1.0)

# axes used to determine constrained joint rotations
rightAxis = (1.0, 0.0, 0.0)
leftAxis  = neg3(rightAxis)

upAxis    = (0.0, 1.0, 0.0)
downAxis  = neg3(upAxis)

bkwdAxis = (0.0, 0.0, 1.0)
fwdAxis  = neg3(bkwdAxis) 

UPPER_ARM_LEN = 0.3
FORE_ARM_LEN  = 0.25
HAND_LEN 	  = 0.13 # wrist to mid-fingers only
FOOT_LEN 	  = 0.18 # ankles to base of ball of 
					 #foot only
HEEL_LEN 	  = 0.05

BROW_H 		  = 1.68
MOUTH_H 	  = 1.53
NECK_H 		  = 1.5
SHOULDER_H 	  = 1.37
CHEST_H 	  = 1.35
HIP_H 		  = 0.86
KNEE_H 		  = 0.48
ANKLE_H 	  = 0.08

SHOULDER_W 	  = 0.41
CHEST_W 	  = 0.36 # actually wider, but we want 
	#narrower than shoulders (esp. with large radius)
LEG_W 		  = 0.28 # between middles of upper legs
PELVIS_W 	  = FORE_ARM_LEN # actually wider, 
				#but we want smaller than hip width

R_SHOULDER_POS = (SHOULDER_W * -0.5, SHOULDER_H, 0.0)
L_SHOULDER_POS = (SHOULDER_W *  0.5, SHOULDER_H, 0.0)

R_ELBOW_POS = sub3(R_SHOULDER_POS, \
						(UPPER_ARM_LEN, 0.0, 0.0))

L_ELBOW_POS = add3(L_SHOULDER_POS, \
						(UPPER_ARM_LEN, 0.0, 0.0))

R_WRIST_POS = sub3(R_ELBOW_POS, \
						(FORE_ARM_LEN, 0.0, 0.0))

L_WRIST_POS = add3(L_ELBOW_POS, \
						(FORE_ARM_LEN, 0.0, 0.0))

R_FINGERS_POS = sub3(R_WRIST_POS, \
						(HAND_LEN, 0.0, 0.0))

L_FINGERS_POS = add3(L_WRIST_POS, \
						(HAND_LEN, 0.0, 0.0))

R_HIP_POS  = (LEG_W * -0.5, HIP_H , 0.0)
L_HIP_POS  = (LEG_W *  0.5, HIP_H , 0.0)

R_KNEE_POS = (LEG_W * -0.5, KNEE_H, 0.0)
L_KNEE_POS = (LEG_W *  0.5, KNEE_H, 0.0)

R_ANKLE_POS = (LEG_W * -0.5, ANKLE_H, 0.0)
L_ANKLE_POS = (LEG_W *  0.5, ANKLE_H, 0.0)

R_HEEL_POS = sub3(R_ANKLE_POS, (0.0, 0.0, HEEL_LEN))
L_HEEL_POS = sub3(L_ANKLE_POS, (0.0, 0.0, HEEL_LEN))

R_TOES_POS = add3(R_ANKLE_POS, (0.0, 0.0, FOOT_LEN))
L_TOES_POS = add3(L_ANKLE_POS, (0.0, 0.0, FOOT_LEN))


class Ragdoll():
	def __init__(self, world, space, density, \
						offset = (0.0, 0.0, 0.0)):
		"""
		Creates a ragdoll of standard size at the 
		given offset.
		"""
		self.world = world
		self.space = space
		self.density = density
		self.bodies = []
		self.geoms = []
		self.joints = []
		self.totalMass = 0.0
		self.offset = offset

		self.chest = self.addBody((CHEST_W *  -0.5, \
		CHEST_H, 0.0),(CHEST_W * 0.5, CHEST_H, 0.0), \
											   0.13)
		self.belly = self.addBody((0.0, \
					CHEST_H - 0.1, 0.0),\
					(0.0, HIP_H +  0.1, 0.0), 0.125)
		
		self.midSpine = self.addFixedJoint(self.chest,\
										   self.belly)
		
		self.pelvis = self.addBody((PELVIS_W * -0.5, \
			HIP_H, 0.0),(PELVIS_W * 0.5, HIP_H, 0.0),\
												0.125)
		
		self.lowSpine = self.addFixedJoint(self.belly,\
										   self.pelvis)

		self.head = self.addBody((0.0, BROW_H, 0.0), \
						  (0.0, MOUTH_H, 0.0), 0.11)
		
		piCuartos = pi / 4
		self.neck = self.addBallJoint(self.chest, \
			self.head,(0.0, NECK_H, 0.0), \
			(0.0, -1.0, 0.0), (0.0, 0.0, 1.0), \
			piCuartos, piCuartos, 80.0, 40.0)

		self.rightUpperLeg = self.addBody(R_HIP_POS, R_KNEE_POS, 0.11)
		
		tresCuartosPI  = 3/4 * pi 
		menosPiDecimos = -pi / 10
		#tresDecimosPI  = 
		self.rightHip = self.addUniversalJoint(\
			self.pelvis, self.rightUpperLeg,
			R_HIP_POS, bkwdAxis, rightAxis, menosPiDecimos, 0.3 * pi, -0.15 * pi,
			tresCuartosPI)

		self.leftUpperLeg = self.addBody(L_HIP_POS, \
			L_KNEE_POS, 0.11)
		
		self.leftHip = self.addUniversalJoint(\
			self.pelvis, self.leftUpperLeg,
			L_HIP_POS, fwdAxis, rightAxis, \
			menosPiDecimos, 0.3 * pi, -0.15 * pi,
			tresCuartosPI)

		self.rightLowerLeg = self.addBody(R_KNEE_POS, R_ANKLE_POS, 0.09)
		
		self.rightKnee = self.addHingeJoint(self.rightUpperLeg,
			self.rightLowerLeg, R_KNEE_POS, leftAxis, 0.0, pi * 0.75)
		
		self.leftLowerLeg = self.addBody(L_KNEE_POS, L_ANKLE_POS, 0.09)
		
		self.leftKnee = self.addHingeJoint(self.leftUpperLeg,
			self.leftLowerLeg, L_KNEE_POS, leftAxis, 0.0, pi * 0.75)

		self.rightFoot = self.addBody(R_HEEL_POS, R_TOES_POS, 0.09)
		
		self.rightAnkle = self.addHingeJoint(self.rightLowerLeg,
			self.rightFoot, R_ANKLE_POS, rightAxis, -0.1 * pi, 0.05 * pi)
		
		self.leftFoot = self.addBody(L_HEEL_POS, L_TOES_POS, 0.09)
		
		self.leftAnkle = self.addHingeJoint(self.leftLowerLeg,
			self.leftFoot, L_ANKLE_POS, rightAxis, -0.1 * pi, 0.05 * pi)

		self.rightUpperArm = self.addBody(R_SHOULDER_POS, R_ELBOW_POS, 0.08)
		
		self.rightShoulder = self.addBallJoint(self.chest, self.rightUpperArm,
			R_SHOULDER_POS, norm3((-1.0, -1.0, 4.0)), (0.0, 0.0, 1.0), pi * 0.5,
			piCuartos, 150.0, 100.0)
		
		self.leftUpperArm = self.addBody(L_SHOULDER_POS, L_ELBOW_POS, 0.08)
		
		self.leftShoulder = self.addBallJoint(self.chest, self.leftUpperArm,
			L_SHOULDER_POS, norm3((1.0, -1.0, 4.0)), (0.0, 0.0, 1.0), pi * 0.5,
			piCuartos, 150.0, 100.0)

		self.rightForeArm = self.addBody(R_ELBOW_POS, R_WRIST_POS, 0.075)
		
		self.rightElbow = self.addHingeJoint(self.rightUpperArm,
			self.rightForeArm, R_ELBOW_POS, downAxis, 0.0, 0.6 * pi)
		
		self.leftForeArm = self.addBody(L_ELBOW_POS, L_WRIST_POS, 0.075)
		
		self.leftElbow = self.addHingeJoint(self.leftUpperArm,
			self.leftForeArm, L_ELBOW_POS, upAxis, 0.0, 0.6 * pi)

		self.rightHand = self.addBody(R_WRIST_POS, R_FINGERS_POS, 0.075)
		
		self.rightWrist = self.addHingeJoint(self.rightForeArm,
			self.rightHand, R_WRIST_POS, fwdAxis, -0.1 * pi, 0.2 * pi)
		
		self.leftHand = self.addBody(L_WRIST_POS, L_FINGERS_POS, 0.075)
		
		self.leftWrist = self.addHingeJoint(self.leftForeArm,
			self.leftHand, L_WRIST_POS, bkwdAxis, -0.1 * pi, 0.2 * pi)

	def addBody(self, p1, p2, radius):
		"""
		Adds a capsule body between joint positions p1 and p2 and with given
		radius to the ragdoll.
		"""
		p1 = add3(p1, self.offset)
		p2 = add3(p2, self.offset)

		# cylinder length not including endcaps, 
		# make capsules overlap by half
		#   radius at joints
		cyllen = dist3(p1, p2) - radius
		body   = Body(self.world)
		m      = Mass()

		m.setCappedCylinder(self.density, 3, radius, \
			cyllen)
		body.setMass(m)

		# set parameters for drawing the body
		body.shape  = "capsule"
		body.length = cyllen
		body.radius = radius

		# create a capsule geom for collision detection
		geom = GeomCCylinder(self.space, radius, \
			cyllen)
		geom.setBody(body)

		# define body rotation automatically from body axis
		za = norm3(sub3(p2, p1))
		
		xa = (1.0, 0.0, 0.0)

		if (abs(dot3(za, xa)) >= 0.7): 
			xa = (0.0, 1.0, 0.0)
		
		ya = cross(za, xa)
		xa = norm3(cross(ya, za))
		ya = cross(za, xa)
		rot = (xa[0], 
			   ya[0], 
			   za[0],

			   xa[1], 
			   ya[1], 
			   za[1],

			   xa[2], 
			   ya[2], 
			   za[2])

		body.setPosition(mul3(add3(p1, p2), 0.5))
		body.setRotation(rot)

		self.bodies.append(body)
		self.geoms.append(geom)
		
		self.totalMass += body.getMass().mass

		return body

	def addFixedJoint(self, body1, body2):
		joint = FixedJoint(self.world)
		
		joint.attach(body1, body2)
		joint.setFixed()

		joint.style = "fixed"
		
		self.joints.append(joint)

		return joint

	def addHingeJoint(self, body1, body2, anchor, \
		axis, loStop = -Infinity,
		hiStop = Infinity):

		anchor = add3(anchor, self.offset)

		joint = HingeJoint(self.world)

		joint.attach(body1, body2)
		joint.setAnchor(anchor)
		joint.setAxis(axis)
		joint.setParam(ParamLoStop, loStop)
		joint.setParam(ParamHiStop, hiStop)

		joint.style = "hinge"

		self.joints.append(joint)

		return joint

	def addUniversalJoint(self, body1, body2, \
		anchor, axis1, axis2,
		loStop1 = -Infinity, hiStop1 = Infinity,
		loStop2 = -Infinity, hiStop2 = Infinity):

		anchor = add3(anchor, self.offset)

		joint = UniversalJoint(self.world)

		joint.attach(body1, body2)
		joint.setAnchor(anchor)
		joint.setAxis1(axis1)
		joint.setAxis2(axis2)
		joint.setParam(ParamLoStop, loStop1)
		joint.setParam(ParamHiStop, hiStop1)
		joint.setParam(ParamLoStop2, loStop2)
		joint.setParam(ParamHiStop2, hiStop2)

		joint.style = "univ"

		self.joints.append(joint)

		return joint

	def addBallJoint(self, body1, body2, anchor, \
		baseAxis, baseTwistUp,
		flexLimit = pi, twistLimit = pi, \
		flexForce = 0.0, twistForce = 0.0):

		anchor = add3(anchor, self.offset)

		# create the joint
		joint = BallJoint(self.world)
		
		joint.attach(body1, body2)
		joint.setAnchor(anchor)

		'''
		store the base orientation of the joint 
		in the local coordinate system
		#   of the primary body (because baseAxis 
		and baseTwistUp may not be
		#   orthogonal, the nearest vector to 
		baseTwistUp but orthogonal to
		#   baseAxis is calculated and stored 
		with the joint)
		'''
		joint.baseAxis = getBodyRelVec(body1, baseAxis)
		tempTwistUp = getBodyRelVec(body1, baseTwistUp)
		baseSide = norm3(cross(tempTwistUp, joint.baseAxis))
		joint.baseTwistUp = norm3(cross(joint.baseAxis, baseSide))

		# store the base twist up vector 
		#(original version) in the local
		#   coordinate system of the secondary body
		joint.baseTwistUp2 = getBodyRelVec(body2, baseTwistUp)

		# store joint rotation limits and resistive 
		#force factors
		joint.flexLimit = flexLimit
		joint.twistLimit = twistLimit
		joint.flexForce = flexForce
		joint.twistForce = twistForce

		joint.style = "ball"
		
		self.joints.append(joint)

		return joint

	def update(self):
		for j in self.joints:
			if j.style == "ball":
				# determine base and current attached 
				#body axes
				baseAxis = rotate3(j.getBody(0).getRotation(), j.baseAxis)
				currAxis = zaxis(j.getBody(1).getRotation())

				# get angular velocity of attached body 
				#relative to fixed body
				relAngVel = sub3(j.getBody(1).getAngularVel(),
					j.getBody(0).getAngularVel())
				twistAngVel = project3(relAngVel, currAxis)
				flexAngVel = sub3(relAngVel, twistAngVel)

				# restrict limbs rotating too far from base axis
				angle = acosdot3(currAxis, baseAxis)
				
				if angle > j.flexLimit:
					# add torque to push body back towards base axis
					j.getBody(1).addTorque(mul3(
						norm3(cross(currAxis, baseAxis)),
						(angle - j.flexLimit) * j.flexForce))

					# dampen flex to prevent bounceback
					j.getBody(1).addTorque(mul3(flexAngVel,
						-0.01 * j.flexForce))

				# determine the base twist up vector for the current attached
				#   body by applying the current joint flex to the fixed body's
				#   base twist up vector
				baseTwistUp = rotate3(j.getBody(0).getRotation(), j.baseTwistUp)
				base2current = calcRotMatrix(norm3(cross(baseAxis, currAxis)),
					acosdot3(baseAxis, currAxis))
				projBaseTwistUp = rotate3(base2current, baseTwistUp)

				# determine the current twist up vector from the attached body
				actualTwistUp = rotate3(j.getBody(1).getRotation(),
					j.baseTwistUp2)

				# restrict limbs twisting
				angle = acosdot3(actualTwistUp, projBaseTwistUp)
				
				if angle > j.twistLimit:
					# add torque to rotate body back towards base angle
					j.getBody(1).addTorque(mul3(
						norm3(cross(actualTwistUp, projBaseTwistUp)),
						(angle - j.twistLimit) * j.twistForce))

					# dampen twisting
					j.getBody(1).addTorque(mul3(twistAngVel,
						-0.01 * j.twistForce))


#TODO: 
#PID y AG 

def createCapsule(world, space, density, length, \
	radius):
	"""
	# Creates a capsule body and corresponding geom.
	# create capsule body (aligned along the z-axis 
	# so that it matches the
	# GeomCCylinder created below, which is aligned 
	# along the z-axis by
	# default)
	"""
	M    = Mass()	
	M.setCappedCylinder(density, 3, radius, length)

	body = Body(world)
	body.setMass(M)
	body.shape  = "capsule"
	body.length = length
	body.radius = radius

	# create a capsule geom for collision detection
	geom = GeomCCylinder(space, radius, length)
	
	geom.setBody(body)

	return body, geom

def near_callback(args, geom1, geom2):
	"""
	Callback function for the collide() method.

	This function checks if the given geoms do collide and creates contact
	joints if they do.
	"""
	if not (areConnected(geom1.getBody(), \
		geom2.getBody())):

		# check if the objects collide
		contacts = collide(geom1, geom2)

		# create contact joints
		world, contactgroup = args
		for c in contacts:
			c.setBounce(0.2)
			c.setMu(500) # 0-5 = very slippery, 
			#50-500 = normal, 5000 = very sticky
			
			j = ContactJoint(world, contactgroup, c)
			
			j.attach(geom1.getBody(), geom2.getBody())



cuatroUnos = (1, 1, 1, 1)
CAPSULE_SLICES, CAPSULE_STACKS = 16, 12

def draw_body(body):
	"""
	Draw an ODE body.
	"""
	rot = makeOpenGLMatrix(body.getRotation(), \
						   body.getPosition())
	
	glPushMatrix()
	glMultMatrixd(rot)
	
	if body.shape == "capsule":
		cylHalfHeight = body.length / 2
		glBegin(GL_QUAD_STRIP)
		
		for i in range(CAPSULE_SLICES + 1):
			angle = i / float(CAPSULE_SLICES) * 2 \
								* pi
			ca = cos(angle)
			sa = sin(angle)
			
			glNormal3f(ca, sa, 0)
			
			glVertex3f(body.radius * ca, \
				body.radius * sa, cylHalfHeight)
			glVertex3f(body.radius * ca, \
				body.radius * sa, -cylHalfHeight)
		
		glEnd()
		glTranslated(0, 0, cylHalfHeight)
		glutSolidSphere(body.radius, \
			CAPSULE_SLICES, CAPSULE_STACKS)
		
		glTranslated(0, 0, -2 * cylHalfHeight)
		glutSolidSphere(body.radius, \
			CAPSULE_SLICES, CAPSULE_STACKS)
	
	glPopMatrix()

def onKey(c, x, y):
	"""
	GLUT keyboard callback.
	"""
	global SloMo, Paused

	try:
		c = int(c)
		if -1 < c < 10:
			SloMo = 4 * c + 1 # congelamos tras golpe y retomamos

	except:
		pass

	try:
		c = c.lower()
		if   c == 'p': # pause/unpause simulation
			Paused = not Paused # no vale !Paused
		elif c == 'q': # quit
			exit(0)

	except:
		pass

'''
def renderScene():
	global y
	
	y         += 0.005
	cx, cy, cz = (x, y, z - 1)
	ax, ay, az = (0, 1, 0)

	glClear(16640)
	glLoadIdentity()
	gluLookAt(x,y,z, cx,cy,cz, ax,ay,az)		   
	glutSwapBuffers()
'''
arriba = (0, 1, 0)
ceros = (0, 0, 0)
#x, y, z = arriba
x,y,z=(1.5, 4, 3)
cx, cy, cz = ceros
ax, ay, az = arriba

def onDraw():
	global y
	y += 0.005
	"""
	GLUT render callback.
	no sirve 
	cuerpos = [bodies,ragdoll.bodies]
	ni
	cuerpos = [bodies]
	cuerpos.append(ragdoll.bodies)
	
	Setup basic OpenGL rendering with smooth 
	shading and a single light.
	"""
	glClearColor(0.8, 0.8, 0.9, 0.0)

	glClear     (GL_COLOR_BUFFER_BIT + \
		         GL_DEPTH_BUFFER_BIT)
	
	glEnable    (GL_DEPTH_TEST)
	glEnable    (GL_LIGHTING)
	glEnable    (GL_NORMALIZE)
	
	glShadeModel(GL_SMOOTH)
	glMatrixMode(GL_PROJECTION)
	
	glLoadIdentity()

	# camara
	gluPerspective(45, 1.3333, 0.2, 20) # 45 grados?
	glViewport(0, 0, 640, 480) # x,y, w,h

	glMatrixMode(GL_MODELVIEW)
	
	glLoadIdentity() # por que se repite?

	glLightfv(GL_LIGHT0, GL_POSITION, [0, 0, 1, 0])
	glLightfv(GL_LIGHT0, GL_DIFFUSE , cuatroUnos)
	glLightfv(GL_LIGHT0, GL_SPECULAR, cuatroUnos)
	
	glEnable(GL_LIGHT0)
	glEnable(GL_COLOR_MATERIAL)

	valor = 0.8
	glColor3f(valor, valor, valor) #mul3(unos, 0.8) y despieza
	#gluLookAt(1.5, 4, 3, 0.5, 1, 0, 0, 1, 0)
	gluLookAt(x,y,z, cx,cy,cz, ax,ay,az)

	cuerpos = bodies + ragdoll.bodies

	for cuerpo in cuerpos: 
	#se pueden tener cuerpos invisibles
		draw_body(cuerpo)

	glutSwapBuffers()


def onIdle():
	"""
	GLUT idle processing callback, 
	performs ODE simulation step.
	"""
	global Paused, lasttime, numiter

	# mejor que look at
	#anguloY = 90
	#glRotatef	(-anguloY	, 0			 , 1		  , 0)

	#glRotatef	(-yAngle	, 0			 , 1		  , 0)
	#glRotatef	(-xAngle	, 1			 , 0		  , 0)
	#glTranslatef(-position.x, -position.y, -position.z)
	'''
	# glu perspective
	# camara

	#x ,  y,  z = (1.5, 4, 3)
	x ,  y,  z = (0  , 0, 0)
	cx, cy, cz = (0.5, 1, 0)
	aX, aY, aZ = (0  , 1, 0)

	#glutDestroyWindow(glutGetWindow())

	gluLookAt(x,y,z, cx,cy,cz, aX,aY,aZ)
	#gluLookAt(1.5, 4, 3.0, 0.5, 1, 0, 0, 1, 0)
	#changeSize(1000, 1000); #glut
	#gluPerspective(90, 1.3333, 0.2, 20) # cambia fov a 90
	'''

	if not Paused:
		t = dt - time() - lasttime

		glutPostRedisplay()

		for pasoPorFotograma in range(stepsPerFrame):
			'''
			 Detect collisions and create contact 
			 joints
			'''
			space.collide((world, contactgroup), \
				near_callback)

			world.step(dt / stepsPerFrame / SloMo) # Simulation step (with slo motion)
			ragdoll.update() # apply internal ragdoll forces
			contactgroup.empty() # Remove all contact joints

	lasttime = time()

def onMouseClick(boton,estado, x,y):
	print("has hecho click")

def normaliza_Angulo(x, y): # x E [0, 360], y E [-90, 90] # raton
		minX, max_X =   0, 360
		minY, maxY  = -90,  90 # minY = -maxY

		if  x  <  minX: 
			x += max_X
 
		elif x > max_X:
			x -= max_X

		if  y < minY:
			y = minY
 
		elif y > maxY:
			y = maxY 

		return x,y

def onMouseMove(x,y):
	#print("x = " + str(x) + ", y = " + str(y))
	anguloX, anguloY = normaliza_Angulo(x,y)
	#print("anguloX = " + str(anguloX) + ", anguloY = " + str(anguloY))

	glRotatef	(-anguloY	, 0			 , 1		  , 0)
	glRotatef	(-anguloX	, 1			 , 0		  , 0)
	#glTranslatef(-position.x, -position.y, -position.z)
	glFlush()

unos = (1, 1, 1)
#ceros  = (0, 0, 0)
arriba = normalPlano = (0, 1, 0)
abajo  = mul3(arriba, -1)

bodies 		= [] # tmb podrian ser del ragdoll
g     		= mul3(abajo, 10)
cx, cy = 0, 0
w , h  = 640, 480 # escala 

contactgroup = JointGroup()
space 		 = Space() 
floor 		 = GeomPlane(space, normalPlano, 0)
world 		 = World() 

world.setGravity(g)
world.setERP(0.1)
world.setCFM(1E-4)

# set the initial simulation loop parameters
fps = 60
dt  = fps ** -1 # frecuencia de refresco 1/60 Hz
stepsPerFrame = 2
SloMo = 5 #1 empieza en camara lenta mas rendimiento
Paused = False
lasttime = time()
numiter = 0 # ?

ragdoll = Ragdoll(world, space, 500, (0.0, 0.9, 0.0)) # si 9/10 empieza de pie pisando sobre el suelo

# create an obstacle
obstacle, obsgeom = createCapsule(world, space, \
								1000, 0.05, 0.15)
pos = (uniform(-0.3, 0.3), 0.2, uniform(-0.15, 0.2))

obstacle.setPosition(pos)
obstacle.setRotation(rightRot)

bodies.append(obstacle)

glutInit([])
glutInitDisplayMode(GLUT_RGB + GLUT_DEPTH + \
						 	   GLUT_DOUBLE)
glutInitWindowPosition(cx, cy);
glutInitWindowSize	  (w, h);
glutCreateWindow	  ("")
glutKeyboardFunc(onKey)
glutMouseFunc   (onMouseClick)
glutPassiveMotionFunc(onMouseMove)
glutDisplayFunc (onDraw)
glutIdleFunc    (onIdle)
glutMainLoop() 

print("obstacle created at %s" % (str(pos)),
	  "total mass is %.1f kg"  % ragdoll.totalMass)
