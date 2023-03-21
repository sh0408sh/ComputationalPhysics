# Pulse propagation simulation by solving a wave equation. Represented by movie

import numpy as np
import matplotlib.pyplot as plt

# Input
xmax=float(input("xmax="))
dx=float(input("dx="))
dt=float(input("dt="))
smax=int(input("smax="))        # max. num. steps of iteration
dsav=int(input("showEvery="))   # show figures every this step 

# Simulation parameters
al=dt/dx   # alpha parameter
a=al**2  
b=2*(1-a)

# Mesh
x=np.arange(0,xmax+dx,dx)      

# Initialize wave variables u and uo
# u is current values and uo is one step older values.
uo=np.exp(-(x-0.2)**2/0.02**2)   # initial pulse at t=0. 
u=np.exp(-(x-0.2-dt)**2/0.02**2) # pulse after one time step advance
ubuf=0*x                         # buffer array for temporary save of var. 
N=len(u)-1                       # last index of the array
u[0]=uo[0]=u[N]=uo[N]=0          # zero boundary condition 

# iteration
s=2    # start from step 2. step 0 and 1 are used for initial conditions
while s<=smax:
	ubuf[1:-1] = u[1:-1]  # save u to buffer
	u[1:-1] = b*u[1:-1] +a*(u[2:]+u[0:-2]) - uo[1:-1] # update u		
	uo[1:-1] = ubuf[1:-1]  # update uo by old u (saved in ubuf) 
      # Shot the figure every 'dsav'th step
	if s%dsav==0:  
		plt.ylim(-1.2,1.2)  # set the ylimit of sub-panels 
		plt.yticks(np.arange(-1.2,1.4,0.4)) # yticks
		plt.plot(x,u); plt.draw(); plt.pause(0.01); plt.clf()
	s+=1  # time step advance

# Show the last shot
plt.show()  

