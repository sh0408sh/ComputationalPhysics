import numpy as np
import random
import matplotlib.pyplot as plt
import plasma2
from mpl_toolkits.mplot3d import Axes3D



#********************* Input Parameters in SI unit ********************* 
# Simulation Domain
nStep = 500
Dt = 0.2
Xlo=-2;  Xup= 2
Ylo=-2;  Yup= 2
Zlo= 0;  Zup= 2


def E(_x, _y, _z):
	#_ex,_ey,_ez = plasma2.ElineXYZ(_x,_y,_z,0,0,0,2)
	#return 0.05*_ex, 0.05*_ey, 0.05*_ez
	return 0.0, 0.0, 0.0

def B(_x, _y, _z):
	_bx,_by,_bz = plasma2.BwireXYZ(_x,_y,_z,0,0,0,2)
	return _bx, _by, _bz


# field line control
ds =0.1
Smax = 60 

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
X[0] = 1.0 
Y[0] = 0.0 
Z[0] = 0.5
Ux[0] = 0.1
Uy[0] = 0.05
#Uy[0] = 0.0
Uz[0] = 0.00

sol1 = plasma2.fLineXYZ(1.0, 0.0, 0.0, ds, Smax, B)
sol2 = plasma2.fLineXYZ(1.0, 0.0, 0.5, ds, Smax, B)
sol3 = plasma2.fLineXYZ(1.0, 0.0, 1.0, ds, Smax, B)

#sol4 = plasma2.fLineXYZ(0.2, 0.0, 0.5, ds, Smax, E)
#sol5 = plasma2.fLineXYZ(0.0, 0.2, 0.5, ds, Smax, E)
#sol6 = plasma2.fLineXYZ(-0.2, 0.0, 0.5, ds, Smax, E)
#sol7 = plasma2.fLineXYZ(0.0, -0.2, 0.5, ds, Smax, E)

fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')
ax.set_xlim(Xlo, Xup)
ax.set_ylim(Ylo, Yup)
ax.set_zlim(Zlo,Zup)
ax.plot3D(sol1[:,0],sol1[:,1],sol1[:,2])
ax.plot3D(sol2[:,0],sol2[:,1],sol2[:,2])
ax.plot3D(sol3[:,0],sol3[:,1],sol3[:,2])
#ax.plot3D(sol4[:,0],sol4[:,1],sol4[:,2])
#ax.plot3D(sol5[:,0],sol5[:,1],sol5[:,2])
#ax.plot3D(sol6[:,0],sol6[:,1],sol6[:,2])
#ax.plot3D(sol7[:,0],sol7[:,1],sol7[:,2])
plt.show()

fig = plt.figure()

while stp<nStep:

	plasma2.moveAF(1,1.0,1.0,Dt,X,Y,Z,Gm,U,Ux,Uy,Uz,E,B)
	Xtr.append(X[0])
	Ytr.append(Y[0])
	Ztr.append(Z[0])

	stp+=1;  T+=Dt
	if stp%10==0:
		ax = fig.add_subplot(111,projection='3d')
		ax.set_xlim(Xlo,Xup); ax.set_ylim(Ylo,Yup); ax.set_zlim(Zlo,Zup)
		ax.plot3D(sol1[:,0],sol1[:,1],sol1[:,2])
		ax.plot3D(sol2[:,0],sol2[:,1],sol2[:,2])
		ax.plot3D(sol3[:,0],sol3[:,1],sol3[:,2])
		#ax.plot3D(sol4[:,0],sol4[:,1],sol4[:,2])
		#ax.plot3D(sol5[:,0],sol5[:,1],sol5[:,2])
		#ax.plot3D(sol6[:,0],sol6[:,1],sol6[:,2])
		#ax.plot3D(sol7[:,0],sol7[:,1],sol7[:,2])
		ax.plot3D(Xtr,Ytr,Ztr)
		plt.draw(); plt.pause(0.01);plt.clf()

ax = fig.add_subplot(111,projection='3d')
ax.set_xlim(Xlo,Xup); ax.set_ylim(Ylo,Yup); ax.set_zlim(Zlo, Zup)
ax.plot3D(sol1[:,0],sol1[:,1],sol1[:,2])
ax.plot3D(sol2[:,0],sol2[:,1],sol2[:,2])
ax.plot3D(sol3[:,0],sol3[:,1],sol3[:,2])
ax.plot3D(Xtr,Ytr,Ztr);plt.show()
