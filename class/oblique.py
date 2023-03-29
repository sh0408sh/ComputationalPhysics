import numpy as np
import matplotlib.pyplot as plt
import emitEMwave as ant

xmax=float(8)
ymax=float(8)
dx=float(0.01)
dy=float(0.01)
dt=float(.005)
f=float(2)
D=float(1)
smax=int(1500)

a=dt/dx;  b=dt/dy
w=2.0*np.pi*f

cntr=int(0.5*ymax)
upper=int(0.5*(ymax+D)/dy)
lower=int(0.5*(ymax-D)/dy)
#upper=int(D/dy)
#lower=int(0)

x=np.arange(0,xmax+dx,dx)
y=np.arange(0,ymax+dy,dy)
Nh = int(xmax/dx/2)

X,Y=np.meshgrid(x,y)
Ex=0*X; Ey=0*X; Ez=0*X
Bx=0*X; By=0*X; Bz=0*X

s=0
while s<smax:
	Ey[:,0]= np.exp(-(y-cntr)**2/(0.2*ymax)**2)*np.sin(w*s*dt) # emission
	#Ey[lower:upper,0]= np.sin(w*s*dt) # hole 
	#ant.emitEMwave(s*dt,Ey[lower:upper,0],(dx,dy),'ovwrt','p',0.8,0,f,0,0,2,(1.0,),30,1)
	Bx[:-1,:-1] +=  -b*(Ez[1:,:-1]-Ez[:-1,:-1])
	By[1:-1,:-1] +=  a*(Ez[1:-1,1:]-Ez[1:-1,:-1])
	Bz[:-1,:-1]+= -a*(Ey[:-1,1:]-Ey[:-1,:-1]) + b*(Ex[1:,:-1]-Ex[:-1,:-1])
	Ex[1:-1,:Nh] +=  b*(Bz[1:-1,:Nh]-Bz[:-2,:Nh])
	Ex[1:-1,Nh:-1] +=  0.25*b*(Bz[1:-1,Nh:-1]-Bz[:-2,Nh:-1])
	Ey[:-1,1:Nh] += -a*(Bz[:-1,1:Nh]-Bz[:-1,:Nh-1])
	Ey[:-1,Nh:-1] += -0.25*a*(Bz[:-1,Nh:-1]-Bz[:-1,Nh-1:-2])
	Ez[1:-1,1:-1]+=  a*(By[1:-1,1:-1]-By[1:-1,:-2]) - b*(Bx[1:-1,1:-1]-Bx[:-2,1:-1])
#	Ez[1:-1,1:-1]+=  a*(By[1:-1,1:-1]-By[1:-1,:-2]) - b*(Bx[1:-1,1:-1]-Bx[:-2,1:-1])
	s+=1


cs=plt.imshow(Ey); plt.colorbar(cs); plt.show()

#plt.plot(x,Ey[cntr, :]);plt.show()
