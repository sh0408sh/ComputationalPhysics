# 1D diffusion solver by Crank-Nicolson method

import numpy as np
import matplotlib.pyplot as plt

# Get the parameters
xmax=float(input("xmax=")) # Domain size
dx=float(input("dx="))   # Mesh size  
dt=float(input("dt="))   # Time step 
D=float(input("D="))   # Diffusion coefficient
imax=int(input("imax="))  # max step to run

# Parameters
a=D*dt/dx**2
b=2*(1+a)
g=2*(1-a)
Nx=int(xmax/dx)

# Variables
x=np.arange(0,Nx+1); x = x*dx
d=0.0*x+b   # Initial diagonal element array
f=0.0*x

# Initial u
u=np.sin(np.pi*x/xmax)+0.3*np.sin(2.0*np.pi*x/xmax)+0.2*np.sin(3.0*np.pi*x/xmax) + 0.1*np.sin(4.0*np.pi*x/xmax)

# Boundary condition
u[0]=u[Nx]=0.0

# Show initial u
plt.ylim(0.0,1.3); plt.yticks(np.arange(0.0,1.4,0.2))
plt.plot(x,u); plt.show()

# Triangularize of the matrix - transform of the diagonal element
uu=-1.0*a   # upper triangluar element
ll=-1.0*a   # lower triangular element
j=0
while j < Nx :
	d[j+1] -= uu*(ll/d[j]);  j+=1


# Main loop
i=0
while i<imax:
	f[1:-1]=g*u[1:-1]+a*u[0:-2]+a*u[2:]
     # transform of f-term 
	j=0
	while j<Nx-1:
		f[j+1] -= f[j]*(ll/d[j]); j+=1 
     # Backsubstitution
	j=Nx-1
	while j>0: 
		u[j] = (f[j] - uu*u[j+1])/d[j];  j -=1
	if i%10==0:
		plt.ylim(0.0,3.3); plt.yticks(np.arange(0.0,3.3,1.1))
		plt.plot(x,u); plt.draw(); plt.pause(0.01); plt.clf()
	i+=1

plt.show()
