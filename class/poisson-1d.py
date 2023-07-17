import numpy as np
import matplotlib.pyplot as plt

# Get some parameters
Nx=int(input("Num. mesh in x="))  # num. meshes
imax=int(input("imax="))   # max step to iterate
tol=float(input("tolerance="))   # error tolerence. something like 1e-5

# Setup arrays. 
x=np.arange(0,Nx+1)  # x-array
u=0.0*x              # u-array, the potential
resid=0.0*x

# Jacobi spectral radius
rjac=np.cos(np.pi/Nx)

# rectangular function for the charge density. Used two Heaviside step func.
s=np.heaviside(x/Nx-0.4,1)*np.heaviside(0.6-x/Nx,1)
plt.plot(x,s);  plt.show()


# initial error 
errf=sum(np.fabs(s))
err= errf  # initial error

# Coefficients
a=1.0
b=1.0
c=-2.0

# initial overrelaxation factor 
omega=1

i=0
while err > (tol*errf) and i<imax: 
      # even calculation
	resid[1:-1:2] = a*u[2::2] +b*u[0:-2:2] + c*u[1:-1:2] - s[1:-1:2]
	u[1:-1:2] -= omega*resid[1:-1:2]/c
	err=sum(np.fabs(resid))
	omega=1/(1-0.5*rjac**2) if i==1 else 1/(1-0.25*rjac**2*omega)

      # odd calculation
	resid[2:-1:2] = a*u[3::2]+b*u[1:-2:2]+c*u[2:-1:2] - s[2:-1:2]
	u[2:-1:2] -= omega*resid[2:-1:2]/c
	err +=sum(np.fabs(resid))
	omega = 1/(1-0.25*rjac**2*omega)

	i+=1

print("Number of iteration=",i)
print("Relative Error=",err/errf)
plt.plot(x,u);  plt.show()
