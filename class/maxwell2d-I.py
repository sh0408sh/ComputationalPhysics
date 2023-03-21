import numpy as np
import matplotlib.pyplot as plt

xmax=float(input("xmax="))
ymax=float(input("ymax="))
dx=float(input("dx="))
dy=float(input("dy="))
dt=float(input("dt="))
f=float(input("f="))
D=float(input("D="))
smax=int(input("smax="))

a=dt/dx;  b=dt/dy
w=2.0*np.pi*f

cntr=0.5*ymax
upper=int(0.5*(ymax+D)/dy)
lower=int(0.5*(ymax-D)/dy)

x=np.arange(0,xmax+dx,dx)
y=np.arange(0,ymax+dy,dy)

X,Y=np.meshgrid(x,y)
Ex=0*X; Ey=0*X; Ez=0*X
Bx=0*X; By=0*X; Bz=0*X

s=0
while s<smax:
	#Ey[:,0]= np.exp(-(y-cntr)**2/(0.2*ymax)**2)*np.sin(w*s*dt) # emission 
	Ey[lower:upper,0]= np.sin(w*s*dt) # hole 
	Bx[:-1,:-1] +=  -b*(Ez[1:,:-1]-Ez[:-1,:-1])
	By[1:-1,:-1] +=  a*(Ez[1:-1,1:]-Ez[1:-1,:-1])
	Bz[:-1,:-1]+= -a*(Ey[:-1,1:]-Ey[:-1,:-1]) + b*(Ex[1:,:-1]-Ex[:-1,:-1])
	Ex[1:-1,:-1] +=  b*(Bz[1:-1,:-1]-Bz[:-2,:-1])
	Ey[:-1,1:-1] += -a*(Bz[:-1,1:-1]-Bz[:-1,:-2])
	Ez[1:-1,1:-1]+=  a*(By[1:-1,1:-1]-By[1:-1,:-2]) - b*(Bx[1:-1,1:-1]-Bx[:-2,1:-1])
	s+=1

cs=plt.imshow(Ey); plt.colorbar(cs); plt.show()
