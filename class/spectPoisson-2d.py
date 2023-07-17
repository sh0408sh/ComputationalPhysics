import numpy as np
from scipy.fftpack import fft, ifft
from scipy.fftpack import fft2, ifft2
import matplotlib.pyplot as plt

xmax=float(input("xmax="))
Nx=int(input("Num. x meshes="))  # power of 2
ymax=float(input("ymax="))
Ny=int(input("Num. y meshes="))  # power of 2
Nxhf=int(Nx/2)
Nyhf=int(Ny/2)

# Create domains
dx=xmax/(1.0*Nx)
dy=ymax/(1.0*Ny)
x= np.arange(0.0,Nx);  x*=dx
y= np.arange(0.0,Ny);  y*=dy

# Create kx, ky arrays
nx= np.arange(0,Nx)
ny= np.arange(0,Ny)
kx= np.zeros(Nx)
ky= np.zeros(Ny)
kx[:Nxhf]=2.0*np.pi*nx[:Nxhf]/(1.0*Nx*dx)
kx[Nxhf:]=2.0*np.pi*(nx[Nxhf:]-Nx)/(1.0*Nx*dx)
ky[:Nyhf]=2.0*np.pi*ny[:Nyhf]/(1.0*Ny*dy)
ky[Nyhf:]=2.0*np.pi*(ny[Nyhf:]-Ny)/(1.0*Ny*dy)

# create a meshgrid
X,Y=np.meshgrid(x,y)

# source
Z=np.exp(-((X/xmax-0.5)**2+(Y/ymax-0.5)**2)/0.04)*np.sin(16.0*np.pi*X/xmax)*np.cos(32*np.pi*Y/ymax)
plt.imshow(Z,origin='lower');  plt.show()


Zh=fft2(Z)
for i in nx:
	for j in ny:
		if i>0 or j>0:
			Zh[i,j] /= -(kx[i]**2+ky[j]**2)
Zh[0,0]=0.0

Z=ifft2(Zh)

plt.imshow(Z.real,origin='lower'); plt.show()
