import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integ


h=float(input("step="))   # time step
nStep=int(input("max step="))  # num. steps to run
m1 = float(input("m1="))   # mass of object 1
m2 = float(input("m2="))   # mass of object 2
m3 = float(input("m3="))   # mass of object 3

# A vector function. y is a two-compo vector, i.e. (x,v).
# f returns two-compo vector, i.e. (f1,f2)=(.
def F(z,t):  
	x1,y1,x2,y2,x3,y3,v1x,v1y,v2x,v2y,v3x,v3y=z  

	r12=np.sqrt((x1-x2)**2+(y1-y2)**2)
	r23=np.sqrt((x2-x3)**2+(y2-y3)**2)
	r31=np.sqrt((x3-x1)**2+(y3-y1)**2)

	f1x=-m2*(x1-x2)/r12**3 - m3*(x1-x3)/r31**3
	f1y=-m2*(y1-y2)/r12**3 - m3*(y1-y3)/r31**3

	f2x=-m1*(x2-x1)/r12**3 - m3*(x2-x3)/r23**3
	f2y=-m1*(y2-y1)/r12**3 - m3*(y2-y3)/r23**3

	f3x=-m1*(x3-x1)/r31**3 - m2*(x3-x2)/r23**3
	f3y=-m1*(y3-y1)/r31**3 - m2*(y3-y2)/r23**3

	return [v1x,v1y,v2x,v2y,v3x,v3y,f1x,f1y,f2x,f2y,f3x,f3y]   

t=np.arange(0,h*nStep,h)  # an array of discretized time 
z0=[0.5,0,  -1.5,0, 0,0,  0,4.3, 0,-1.5,  0,0]  # initial conditions: x0=1, v0=0
sol = integ.odeint(F,z0,t) # sol is [[x0,y0,vx0,vy0],..]

#plt.plot(t,sol); plt.show()
plt.plot(sol[:,0],sol[:,1])
plt.plot(sol[:,2],sol[:,3])
plt.plot(sol[:,4],sol[:,5]);plt.show()

plt.plot(sol[:,0]-sol[:,4],  sol[:,1]-sol[:,5]);  
plt.plot(sol[:,2]-sol[:,4],  sol[:,3]-sol[:,5]);  plt.show()

