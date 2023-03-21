import numpy as np
import matplotlib.pyplot as plt

xmax=float(input("xmax="))
dx=float(input("dx="))
dt=float(input("dt="))
f=float(input("f="))
smax=int(input("smax="))

a=dt/dx
w=2.0*np.pi*f

x=np.arange(0,xmax+dx,dx)
c=int(0.5*xmax/dx)

Ey=0*x; Ez=0*x
By=0*x; Bz=0*x

s=0
while s<smax:
	By[:-1]  +=  a*(Ez[1:]-Ez[:-1])
	Bz[:-1]  += -a*(Ey[1:]-Ey[:-1]) 
	Ey[1:-1] += -a*(Bz[1:-1]-Bz[0:-2]) 
	Ey[c] += -dt*np.sin(w*s*dt)
	Ez[1:-1] +=  a*(By[1:-1]-By[0:-2])
	s+=1

plt.plot(x,Ey); plt.show()
