import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integ
from scipy.fftpack import fft,ifft

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

plt.plot(t,sol); plt.show()

plt.plot(sol[:,0],sol[:,1]); plt.show()

#zh=fft(sol[:,2])
#plt.plot(abs(zh)); plt.show()

