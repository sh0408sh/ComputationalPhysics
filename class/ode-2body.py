import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integ


h=float(input("step="))   # time step
nStep=int(input("max step="))  # num. steps to run
m1 = float(input("m1="))   # mass of object 1
m2 = float(input("m2="))   # mass of object 2

# A vector function. y is a two-compo vector, i.e. (x,v).
# f returns two-compo vector, i.e. (f1,f2)=(.
def F(z,t):  
	x1,y1,x2,y2,v1x,v1y,v2x,v2y=z  
	r=np.sqrt((x1-x2)**2+(y1-y2)**2)
	f1x=-m2*(x1-x2)/r**3; f1y=-m2*(y1-y2)/r**3
	f2x=-m1*(x2-x1)/r**3; f2y=-m1*(y2-y1)/r**3
	return [v1x,v1y,v2x,v2y,f1x,f1y,f2x,f2y]   

t=np.arange(0,h*nStep,h)  # an array of discretized time 
z0=[0.5,0,  0,0,  0,4.3, 0,0]  # initial conditions: x0=1, v0=0
sol = integ.odeint(F,z0,t) # sol is [[x0,y0,vx0,vy0],..]

#plt.plot(t,sol); plt.show()
plt.plot(sol[:,0],sol[:,1])
plt.plot(sol[:,2],sol[:,3]); plt.show()
plt.plot(sol[:,0]-sol[:,2],  sol[:,1]-sol[:,3]);  plt.show()

