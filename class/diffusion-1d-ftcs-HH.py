# 1D diffusion solver by FTCS method
# Suggested parameters: xmax=1, dx=0.01, dt=0.001, D=0.05, imax=1000

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
Nx=int(xmax/dx)

# Variables
x=np.arange(0,Nx+1); x = x*dx
u=0.0*x

# initial shape of u
u=np.sin(np.pi*x/xmax)+0.3*np.sin(2.0*np.pi*x/xmax)+0.2*np.sin(3.0*np.pi*x/xmax) + 0.1*np.sin(4.0*np.pi*x/xmax)

# Show the initial shape of u
plt.ylim(0.0,1.3); plt.yticks(np.arange(0,1.4,0.2)); plt.plot(x,u); plt.show()

# Boundary condition (zero)
u[0]=u[Nx]=0.0

# Main loop
i=0
while i<imax:
	u[1:-1] += a*(u[2:]-2*u[1:-1]+u[0:-2])
	if i%10==0:
		plt.ylim(0.0,1.3)
		plt.yticks(np.arange(0.0,1.4,0.2))
		plt.plot(x,u); plt.draw(); plt.pause(0.01); plt.clf()
	i+=1

plt.show()
