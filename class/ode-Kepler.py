import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integ


h=float(input("step="))   # time step
nStep=int(input("max step="))  # num. steps to run

# A vector function. y is a two-compo vector, i.e. (x,v).
# f returns two-compo vector, i.e. (f1,f2)=(.
def F(t,z):  
	x,y,vx,vy=z  
	r=np.sqrt(x**2+y**2)
	fx=-x/r**3; fy=-y/r**3
	return [vx,vy,fx,fy]   

t=np.arange(0,h*nStep,h)  # an array of discretized time 
z0=[1.0, 0.0, 0.00, 0.5]  # initial conditions: x0=1, v0=0
#sol = integ.odeint(F,z0,t) # sol is [[x0,y0,vx0,vy0],..]
sol = integ.solve_ivp(F,[0,10],z0,method='RK45',t_eval=t) # sol is [[x0,y0,vx0,vy0],..]

#plt.plot(t,sol); plt.show()
#plt.plot(sol[:,0],sol[:,1]); plt.show()
plt.plot(sol.y[0],sol.y[1]); plt.show()

