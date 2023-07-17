import random
import numpy as np
import matplotlib.pyplot as plt
from numba import njit

N=300        # number of molecules
smax=1000   # number of iterations 
T=0.1         # temperature (average kinetic energy, 0.5v^2)
dt = 0.002
L = 1.0     # length of the simulation box
sgm=0.2
eps=0.05

##### Make a plot of LJ potential #####
d_ave = L/N**(1/3)  # average distance between particles
x = np.arange(d_ave*0.01, 1.1*d_ave,0.001)
U = 4*sgm*((eps/x)**12-(eps/x)**6)
plt.plot(x,U); plt.plot(d_ave,0,'ro')
plt.plot([0,1.1*d_ave],[-1.5*T,-1.5*T])
plt.xlim(0,1.1*d_ave)
plt.ylim(-1.7*T,0.5*T)
plt.show()
#######################################

R=np.zeros([N,3])    
Ro=np.zeros([N,3])    
Rtmp=np.zeros([N,3])    
V=np.zeros([N,3])    
F = np.zeros([N,3])

@njit
def F_LJ(_N,_R,_F):
	_e = eps**2
	Utot = 0
	for i in range(_N):
		_F[i,] *= 0.0

		for j in range(_N):
			if j != i:
				X = _R[i,] - _R[j,]
				r2 = np.sum(_R[i,]*_R[i,])
				pli = (_e/r2)**6
				vdw = (_e/r2)**3
				Utot += 4*sgm*(pli-vdw)
				fac = 24*sgm*(2.0*pli-vdw**3)/r2
				_F[i,] += fac*X 
	return 0.5*Utot


# initial velocity distributon
for i in range(N): 
#	_vx = 0.5*(random.random()-0.5)
#	_vy = 0.5*(random.random()-0.5)
#	_vz = 0.5*(random.random()-0.5)
	_vx = random.gauss(0, np.sqrt(2*T)) 
	_vy = random.gauss(0, np.sqrt(2*T)) 
	_vz = random.gauss(0, np.sqrt(2*T)) 
	_x = L*(0.5*random.random() + 0.25)
	_y = L*(0.5*random.random() + 0.25)
	_z = L*(0.5*random.random() + 0.25)
	Ro[i,] = (_x,_y,_z)
	V[i,] = (_vx,_vy,_vz)

plt.plot(Ro[:,0],Ro[:,1],'ro',ms=2); plt.show()
plt.plot(V[:,0],V[:,1],'ro',ms=2); plt.show()

# the 1st step
R = Ro + dt*V + 0.5*dt**2*F
V = (R-Ro)/dt

tarr=[]
Earr=[]
for s in range(smax):
	U = F_LJ(N,R,F)
	K=0   # total kinetic energy
	for i in range(N):
		Rtmp[i,] = R[i,]
		R[i,] = 2*R[i,] - Ro[i,] + dt**2*F[i,]
		Ro[i,] = Rtmp[i,]

		for c in range(3):
			if R[i,c] > L:  
				R[i,c] = 2*L - R[i,c]  
				Ro[i,c] = 2*L - Ro[i,c]
			elif R[i,c] < 0:  
				R[i,c] = -R[i,c]  
				Ro[i,c] = -Ro[i,c]
		V[i,] = (R[i,] - Ro[i,])/dt

		ke = 0.5*np.sum(V[i,]*V[i,])	
		K += ke
		
		if np.sqrt(2*ke)*dt > L:  print(np.sqrt(2*ke))
	tarr.append(s)
	Earr.append(K+U)

#	plt.plot(R[:,0],R[:,1],'ro',ms=2)
#	plt.xlim(0,L)
#	plt.ylim(0,L)
	plt.plot(tarr,Earr)
	plt.draw()
	plt.pause(0.1)
	plt.clf()
plt.show()

plt.plot(tarr,Karr); plt.show()
#plt.ioff()
#plt.plot(x[:,0],x[:,1],'ro',ms=1.0)
#plt.xlim(mn,mx)
#plt.ylim(mn,mx)
plt.show()
