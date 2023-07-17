import random
import numpy as np
import matplotlib.pyplot as plt
from numba import njit

N=300        # number of molecules
smax=200   # number of iterations 
T=0.1         # temperature (average kinetic energy, 0.5v^2)
dt = 0.02
L = 1.0     # length of the simulation box
A=-3
B=0.02
S=0.05
eps=0.05

##### Make a plot of potential #####
d_ave = L/N**(1/3)  # average distance between particles
x = np.arange(d_ave*0.01, 1.1*d_ave,0.001)
U = A*(x-B)*np.exp(-x**2/S**2)
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
def F_MH(_N,_R,_F):
	_e = eps**2
	Utot = 0
	for i in range(_N):
		_F[i,] *= 0.0

		for j in range(_N):
			if j != i:
				X = _R[i,] - _R[j,]
				r2 = np.sum(_R[i,]*_R[i,])
				r = np.sqrt(r2)
				ex = np.exp(-r2/S**2)
				Utot += A*(r-B)*ex
				_F[i,] += -A*(1 - 2*(r-B)*r/S**2)*ex*(X/r)
	return 0.5*Utot


# initial velocity distributon
for i in range(N): 
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
	U = F_MH(N,R,F)
	K=0   # total kinetic energy
	for i in range(N):
		V[i,] += dt*F[i,]
		R[i,] += dt*V[i,]

		for c in range(3):
			if R[i,c] > L:  
				R[i,c] = 2*L - R[i,c]  
				V[i,c] *= -1
			elif R[i,c] < 0:  
				R[i,c] = -R[i,c]  
				V[i,c] *= -1

		ke = 0.5*np.sum(V[i,]*V[i,])	
		K += ke
		
		if np.sqrt(2*ke)*dt > L:  print(np.sqrt(2*ke))
	tarr.append(s)
	Earr.append(K+U)

	plt.plot(R[:,0],R[:,1],'ro',ms=2)
	plt.xlim(0,L)
	plt.ylim(0,L)
#	plt.plot(tarr,Earr)
	plt.draw()
	plt.pause(0.1)
	plt.clf()
plt.show()

plt.plot(tarr,Earr); plt.show()
