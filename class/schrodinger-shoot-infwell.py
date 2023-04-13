import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integ


h=0.005
t=np.arange(0,1.0,h)

def f(y,t,E):
	x,v=y
	return [v,-2.0*E*x]

y0=[0.0,1.0]

i=0
while i==0:
	E=float(input("E="))
	sol = integ.odeint(f,y0,t,args=(E,))
	plt.plot(t,sol[:,0]); plt.show()

