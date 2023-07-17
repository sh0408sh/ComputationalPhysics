import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integ
from mpl_toolkits.mplot3d import Axes3D

s=10
r=float(input("r="))
b=2.66667

h=float(input("step="))
nStep=int(input("max step="))

def f(q,t):
	x,y,z=q
	dxdt=s*(y-x)
	dydt=r*x-y-x*z
	dzdt=x*y-b*z
	return [dxdt,dydt,dzdt]


t=np.arange(0,h*nStep,h)

q0=[1.0, 1.0, 1.0]
sol = integ.odeint(f,q0,t)

fig=plt.figure()
ax=fig.add_subplot(111,projection='3d')
ax.set_xlim(-20,20); ax.set_ylim(-30,30); ax.set_zlim(0,50)
ax.plot3D(sol[:,0],sol[:,1],sol[:,2]); plt.show()



