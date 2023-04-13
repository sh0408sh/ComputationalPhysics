import numpy as np
import matplotlib.pyplot as plt
import emitEMwave as ant
import emitEMwave_update as antu
import math as mt
import time as tm
import refraction as rf

#set the variables

xmax=float(4)
ymax=float(10)
dx=float(0.02)
dy=float(0.02)
dt=float(0.01)
f=float(5)
D=float(3)

#smax=int(2*xmax/(dx*np.cos(phi)))
dsav = 50

n1 = np.sqrt(2)
n2 = 1

inc_angle = np.arange(30.0,60.0,1.0)
#R = inc_angle*0
phi = 45.0*np.pi/180.0
smax=int(2.5*xmax/(dx*np.cos(phi)))
R = rf.refraction(xmax, ymax, dx, dy, dt,f,D,smax,phi,1,n1,n2,0.5)
print(R)

"""
phi = 30.0*np.pi/180.0
R = rf.refraction(xmax, ymax, dx, dy, dt,f,D,smax,phi,1,n1,n2,0.5)
print(R)
"""
"""
plt.plot(inc_angle,R)
plt.xlabel('incident angle')
plt.ylabel('Reflectance')
plt.show()
"""
#start = tm.time()
#print(rf.refraction(xmax, ymax, dx, dy, dt,f,D,smax,np.pi/4,1,n1,n2,0.3))
#end = tm.time()

#print(end-start)