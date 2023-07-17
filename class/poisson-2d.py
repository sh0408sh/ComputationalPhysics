import numpy as np
import matplotlib.pyplot as plt

# Get some parameters
Nx=int(input("Num. mesh in x="))  # num. meshes in x
xmax=float(input("xmax="))        # Max x
Ny=int(input("Num. mesh in y="))  # num. meshes in y
ymax=float(input("ymax="))        # Max y
imax=int(input("imax="))   # max step to iterate
tol=float(input("tolerance="))   # error tolerence. something like 1e-5

# Setup arrays and parameters. 
dx=xmax/Nx    # mesh size in x
dy=ymax/Ny    # mesh size in y
x=np.arange(0,Nx+1)  # x-array
y=np.arange(0,Ny+1)  # y-array
X,Y=np.meshgrid(x,y) # Create a 2D mesh grid
u=0.0*X              # u-array, the potential
resid=0.0*X          # residue 

# Jacobi spectral radius
rjac=(np.cos(np.pi/Nx) + (dx/dy)**2*np.cos(np.pi/Ny))/(1+(dx/dy)**2)

# Gaussian source
s=np.exp(-((X/Nx-0.5)**2+(Y/Ny-0.5)**2)/0.04)

# plot of Charge distribution 
extn=[0,xmax,0,ymax]  # axis range
cs=plt.imshow(s,extent=extn); plt.colorbar(cs); plt.show()

# initial error 
errf=sum(sum(np.fabs(s)))
err= errf  # initial error

# Coefficients
a=1.0
b=1.0
c=(dx/dy)**2
d=(dx/dy)**2
e=-2.0-2.0*(dx/dy)**2

# initial overrelaxation factor 
omega=1

i=0; 
while err > (tol*errf) and i<imax: 
      # odd calculation
	j=1; lsw=1
	while j < Nx:
		resid[lsw:-1:2, j] =  a*u[lsw:-1:2, j+1] + b*u[lsw:-1:2, j-1] \
                                    + c*u[lsw+1::2, j] + d*u[lsw-1:-2:2, j] \
                                    + e*u[lsw:-1:2, j] + s[lsw:-1:2, j] 
		u[lsw:-1:2, j] -= omega*resid[lsw:-1:2, j]/e
		err=sum(sum(np.fabs(resid)))
		j+=1
		lsw=3-lsw
	omega=1/(1-0.5*rjac**2) if i==1 else 1/(1-0.25*rjac**2*omega)

      # even calculation
	j=1;  lsw=2
	while j < Nx:
		resid[lsw:-1:2, j] =  a*u[lsw:-1:2, j+1] + b*u[lsw:-1:2, j-1] \
                                    + c*u[lsw+1::2, j] + d*u[lsw-1:-2:2, j] \
                                    + e*u[lsw:-1:2, j] + s[lsw:-1:2, j] 
		u[lsw:-1:2, j] -= omega*resid[lsw:-1:2, j]/e
		err += sum(sum(np.fabs(resid)))
		j+=1
		lsw=3-lsw
	omega = 1/(1-0.25*rjac**2*omega)
	i+=1

# print out run info.
print("Number of iteration=",i)
print("Relative Error=",err/errf)

# plot of the resulting potential
extn=[0,xmax,0,ymax]
cs=plt.imshow(u,extent=extn); plt.colorbar(cs);  plt.show()
