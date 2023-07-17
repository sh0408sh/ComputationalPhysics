import numpy as np
import random
import matplotlib.pyplot as plt
import plasma2



#********************* Input Parameters in SI unit ********************* 
# Simulation Domain
Smax = 400
Dt = 0.1

Ux0 = 0.01
Uy0 = 0.0
Uz0 = 0.0

Xlower = -0.02
Xupper = 0.02
Ylower = -0.05
Yupper = 0.02

def Prf(_x):
	return 1.0

def E(_x, _y, _z):
	return 0.001, 0.0, 0.0

def B(_x, _y, _z):
	return 0.0, 0.0, 1

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

stp=0;  T=0
Ux[0] = Ux0
Uy[0] = Uy0
Uz[0] = Uz0

while stp<Smax:

	plasma2.moveAF(1,1.0,1.0,Dt,X,Y,Z,Gm,U,Ux,Uy,Uz,E,B)
	Xtr.append(X[0])
	Ytr.append(Y[0])

	stp+=1;  T+=Dt
	if stp%10==0:
		plt.ylim(Ylower,Yupper); plt.xlim(Xlower,Xupper)
		plt.plot(Xtr,Ytr)
		plt.draw(); plt.pause(0.01);plt.clf()

plt.plot(Xtr,Ytr); plt.show()
