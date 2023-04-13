import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integ  # to use odeint()

#V=float(input("V="))  # get the dc voltage
#R=float(input("R="))  # get the resistance
#L=float(input("L="))  # get the inductance
V=3
R=1
L=1
a=R/L

h=float(input("step="))   # time step for odeint()
nStep=int(input("max step="))  # num. steps to run odeint()

# define a function f(x) in dx/dt=f(x,t). x corresponds to i in RL circuit
def f(x,t):
	return -a*x+V/L

t=np.arange(0,h*nStep,h)   # prepare an array of discretized time domain
x0=[0.0]   # initial current is zero
x = integ.odeint(f,x0,t)   # call odeint(func, init, time-array). It returns the solution,i.e. current at each time step as an array.  

plt.plot(t,x); plt.show()

