import numpy as np
import random
import matplotlib.pyplot as plt
import plasma3 as plasma2
from mpl_toolkits.mplot3d import Axes3D



#********************* Input Parameters in SI unit ********************* 
# Simulation Domain
nStep = 3000
Dt = 0.3
Xlo = -0.4;  Xup =  0.4
Ylo = -0.4;  Yup =  0.4
Zlo =  0.0;  Zup =  5.0

X0 = 0.2;  Y0 = 0.0;  Z0 = 2.5
Ux0 = 0.02;  Uy0 = 0.0; Uz0 = 0.05

def E(_x, _y, _z):
	return 0.0, 0.0, 0.0

def B(_x, _y, _z):
	_r = np.sqrt(_x**2 + _y**2)
	_br = plasma2.Br(_r,_z,4.0,0.5)+plasma2.Br(_r,_z,1.0,0.5)
	_bx = _br*_x/_r
	_by = _br*_y/_r
	_bz = plasma2.Bz(_r,_z,4.0,0.5)+plasma2.Bz(_r,_z,1.0,0.5)
	return _bx, _by, _bz

# field line control
ds=0.05
Smax = 70
#********************** End of Input parameters *********************** 

X=np.zeros(1)
Y=np.zeros(1)
Z=np.zeros(1)
U=np.zeros(1)
Gm=np.zeros(1)
Ux=np.zeros(1)
Uy=np.zeros(1)
Uz=np.zeros(1)

Xtr=[]
Ytr=[]
Ztr=[]

stp=0;  T=0
X[0] = X0; Y[0] = Y0;  Z[0] = Z0
Ux[0] = Ux0;  Uy[0] = Uy0;  Uz[0] = Uz0

sol1 = plasma2.fLineXYZ(0.1, 0.0, 0.5, ds, Smax, B)
sol2 = plasma2.fLineXYZ(0.0, 0.1, 0.5, ds, Smax, B)
sol3 = plasma2.fLineXYZ(-0.1, 0.0, 0.5, ds, Smax, B)
sol4 = plasma2.fLineXYZ(0.0, -0.1, 0.5, ds, Smax, B)

fig = plt.figure(figsize=(5,10))
ax = fig.add_subplot(111,projection='3d')
#ax.set_box_aspect(aspect=(1,1,2),zoom=1.5)
#ax.view_init(elev=60,azim=-60)
ax.set_xlim(Xlo, Xup)
ax.set_ylim(Ylo, Yup)
ax.set_zlim(Zlo,Zup)
ax.plot3D(sol1[:,0],sol1[:,1],sol1[:,2])
ax.plot3D(sol2[:,0],sol2[:,1],sol2[:,2])
ax.plot3D(sol3[:,0],sol3[:,1],sol3[:,2])
ax.plot3D(sol4[:,0],sol4[:,1],sol4[:,2])
plt.show()

fig = plt.figure(figsize=(5,10))

while stp< nStep:

	plasma2.moveAF(1,1.0,1.0,Dt,X,Y,Z,Gm,U,Ux,Uy,Uz,E,B)
	Xtr.append(X[0])
	Ytr.append(Y[0])
	Ztr.append(Z[0])

	stp+=1;  T+=Dt
	if stp%10==0:
		ax = fig.add_subplot(111,projection='3d')
		#ax.set_xlim(-0.4,0.4); ax.set_ylim(-0.4,0.4); ax.set_zlim(0, 5)
		ax.set_xlim(Xlo, Xup); ax.set_ylim(Ylo, Yup); ax.set_zlim(Zlo,Zup)
#		ax.view_init(elev=30,azim=-70)
		#ax.set_box_aspect(aspect=(1,1,2),zoom=1.5)
#		ax.view_init(elev=60,azim=-60)
		ax.plot3D(sol1[:,0],sol1[:,1],sol1[:,2])
		ax.plot3D(sol2[:,0],sol2[:,1],sol2[:,2])
		ax.plot3D(sol3[:,0],sol3[:,1],sol3[:,2])
		ax.plot3D(sol4[:,0],sol4[:,1],sol4[:,2])
		ax.plot3D(Xtr,Ytr,Ztr)
		plt.draw(); plt.pause(0.01);plt.clf()

ax = fig.add_subplot(111,projection='3d')
ax.set_xlim(-0.4,0.4); ax.set_ylim(-0.4,0.4); ax.set_zlim(0, 5)
#ax.view_init(elev=90,azim=0)
ax.plot3D(sol1[:,0],sol1[:,1],sol1[:,2])
ax.plot3D(sol2[:,0],sol2[:,1],sol2[:,2])
ax.plot3D(sol3[:,0],sol3[:,1],sol3[:,2])
ax.plot3D(sol4[:,0],sol4[:,1],sol4[:,2])
ax.plot3D(Xtr,Ytr,Ztr)
plt.show()
