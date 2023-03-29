# 2D diffusion solver by ADI method
# Suggested parameters are xmax=ymax=1, dx=dy=0.01, dt=0.005, D=0.05, imax=500

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
#from matplotlib import cm

# Get the parameters
xmax=float(input("xmax=")) # Domain size in x
dx=float(input("dx="))   # Mesh size in x  
ymax=float(input("ymax=")) # Domain size in y
dy=float(input("dy="))   # Mesh size in y 
dt=float(input("dt="))   # Time step 
D=float(input("D="))   # Diffusion coefficient
imax=int(input("imax="))  # max step to run

# Parameters
a=D*dt/dx**2  # alpha
b=2*(1+a)     # beta 
g=2*(1-a)     # gamma
Nx=int(xmax/dx)
Ny=int(ymax/dy)

# Variables
x=np.arange(0,Nx+1); x = x*dx
y=np.arange(0,Ny+1); y = y*dy
X,Y=np.meshgrid(x,y)

# Initial diagonal element array
d1 = 0.0*x +b
d2 = 0.0*y +b   

f=0.0*X     # temporary source term in BTCS part

# Initial u
u=0.0*X
u[int(0.4*Ny):int(0.6*Ny),int(0.4*Nx):int(0.6*Nx)] =1.0 # Source at the center

# Prepare 3D (x,y,t) array for animation
nfrm = int(imax/10)  # number of frames to run
zarr=np.zeros((Ny+1,Nx+1,nfrm)) 
zarr[:,:,:] *=0.0

# Boundary condition
u[0,:]=u[Ny,:]=u[:,0]=u[:,Nx]=0.0

# Triangularize of the matrix - transform of the diagonal element
uu=-1.0*a   # upper triangluar element
ll=-1.0*a   # lower triangular element
j=0
while j < Nx :
	d1[j+1] -= uu*(ll/d1[j]);  j+=1
j=0
while j < Ny :
	d2[j+1] -= uu*(ll/d2[j]);  j+=1

# Main loop
i=0
while i<imax:
   # Every 10 steps, save the data to zarr for later use in animation
	if i%10==0:
		zarr[:,:,int(i/10)] = u[:,:]
   ## half-step advance in x
     # Set the source
	u[int(0.4*Ny):int(0.6*Ny),int(0.4*Nx):int(0.6*Nx)] = 1.0
     # Calculate f_{j,l}^n
	f[1:-1,1:-1]=g*u[1:-1,1:-1]+a*u[0:-2,1:-1]+a*u[2:,1:-1]
     # transform of f 
	j=0
	while j<Nx-1:
		f[1:-1,j+1] -= f[1:-1,j]*(ll/d1[j]); j+=1 
     # Backsubstitution for u^{n+1/2}
	j=Nx-1
	while j>0: 
		u[1:-1,j] = (f[1:-1,j] - uu*u[1:-1,j+1])/d1[j];  j -=1

   ## another half-step advance in y
     # Set the source
	u[int(0.4*Ny):int(0.6*Ny),int(0.4*Nx):int(0.6*Nx)] = 1.0
     # Calculate f_{j,l}^{1/2}
	f[1:-1,1:-1]=g*u[1:-1,1:-1]+a*u[1:-1,0:-2]+a*u[1:-1,2:]
     # transform of f 
	l=0
	while l<Ny-1:
		f[l+1,1:-1] -= f[l,1:-1]*(ll/d2[l]); l+=1 
     # Backsubstitution for u^{n+1}
	l=Ny-1
	while l>0: 
		u[l,1:-1] = (f[l,1:-1] - uu*u[l+1,1:-1])/d2[l];  l -=1
	i+=1


# Animate the result
def update_plot(frm_n,zarr_,plot):
	plot[0].remove()
	plot[0]=ax.plot_surface(X,Y,zarr[:,:,frm_n],cmap="jet")

fig=plt.figure()
ax=fig.add_subplot(111,projection='3d')

plot=[ax.plot_surface(X,Y,zarr[:,:,0], color='0.75', rstride=1,cstride=1)]
ax.set_zlim(0,1.1)
ani=animation.FuncAnimation(fig, update_plot, nfrm,fargs=(zarr,plot),interval=1000/10,repeat=False)

plt.show()
