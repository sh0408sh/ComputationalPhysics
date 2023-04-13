import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integ

#w=float(input("w="))   # get the frequency of HO
w=2

h=float(input("step="))   # time step
nStep=int(input("max step="))  # num. steps to run

# A vector function. y is a two-compo vector, i.e. (x,v).
# f returns two-compo vector, i.e. (f1,f2)=(.
def f(y,t):  
	x,v=y  # by this, x=y[0] (position), v=y[1] (velocity)
	return [v,-w*w*x]   # by this, f returns (f1,f2) where f1=v, f2=-w*w*x

t=np.arange(0,h*nStep,h)  # an array of discretized time 
y0=[1.0,0.0]  # initial conditions: x0=1, v0=0
sol = integ.odeint(f,y0,t) # sol is a 2D array.It is like [[x0,v0],[x1,v1],...]

plt.plot(t,sol); plt.show()
plt.plot(sol[:,0],sol[:,1]); plt.show()
