# 2D Poisson solver, mixed SOR and spectral method.
# Spectral in y and SOR in x.

import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft, ifft


# Get some parameters
xmax=float(input("Xmax="))        # max in x range
Nx=int(input("Num. mesh in x="))  # num. meshes in x
ymax=float(input("Ymax="))        # max in y range
Ny=int(input("Num. mesh in y="))  # num. meshes in y, power of 2
imax=int(input("imax="))   # max step to iterate, for SOR iteration
tol=float(input("tolerance="))   # error tolerence. something like 1e-5. It is for SOR part


# Setup arrays and parameters. 
dx=xmax/Nx           # mesh size in x
dy=ymax/Ny           # mesh size in y
x=np.arange(0.0,Nx); x *=dx  # x-array
y=np.arange(0.0,Ny); y *=dy  # y-array
X,Y=np.meshgrid(x,y) # Create mesh for 2D
u=(0.0+0.0j)*X       # u-array, the potential, 2D. Initialized by complex zero
extn=[0,xmax,0,ymax] # range for imshow()


# Jacobi spectral radius
rjac=np.cos(np.pi/Nx)


# Source of high-k sine in y and smooth Gaussian in x
s=np.sin(64.0*np.pi*Y/ymax) * np.exp(-(X/xmax-0.5)**2/0.04-(Y/ymax-0.5)**2/0.01) 
plt.imshow(s,extent=extn); plt.show()

# Coefficients
a=1.0
b=1.0
c=-2.0


Nhf=int(Ny/2)
k= np.zeros(Ny)
n=np.arange(0,Ny)
k[:Nhf]=n[:Nhf]/(dy*Ny)
k[Nhf:]=(n[Nhf:]-Ny)/(dy*Ny)
k = dx*k

shat=(0.0j)*s
for i in range(Nx):
	shat[:,i]=fft(s[:,i])
shat = (dx**2)*shat

# Main loop
omega=1  # initial overrelaxation factor 
for j in range(Ny):
      # initial error 
	errf=sum(np.abs(shat[j,:]))
	err= errf  # initial error
	c=-2-k[j]**2    # diagonal coefficient 
	i=0
	while err > (tol*errf) and i<imax: 
	      # even calculation
		resid = a*u[j,2::2]+b*u[j,0:-2:2] +c*u[j,1:-1:2] -shat[j,1:-1:2]
		u[j,1:-1:2] -= omega*resid/c
		err=sum(np.abs(resid))
		omega=1/(1-0.5*rjac**2) if i==1 else 1/(1-0.25*rjac**2*omega)

      	# odd calculation
		resid = a*u[j,3::2]+b*u[j,1:-2:2]+c*u[j,2:-1:2] -shat[j,2:-1:2]
		u[j,2:-1:2] -= omega*resid/c
		err +=sum(np.abs(resid))
		omega = 1/(1-0.25*rjac**2*omega)

		i+=1

#	print("Number of iteration=",i)
#	print("Relative Error=",err/errf)

uR=0.0*X
for j in range(Nx):
	tmp=ifft(u[:,j])
	uR[:,j]=tmp.real
cs=plt.imshow(uR,extent=extn);  plt.colorbar(cs);  plt.show()
