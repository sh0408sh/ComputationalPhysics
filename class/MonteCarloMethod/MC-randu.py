import numpy as np
import matplotlib.pyplot as plt

a=65539
b=0
c=2**31

N=int(input("How many random numbers? ="))
x=np.arange(0,N)
y=0.0*x

z=1
for var in x:
	y[var]=z
	z=(a*z+b)%c
	
y=(1/c)*y            # Normalization of y 
plt.plot(x,y,'ro',ms=1); plt.show()	

# Plot the distribution
f=np.arange(0,101)   # distribution function f. Num. bins is 100
f=0.0*f     # initialization of f by zero
for var in y:
	scaly=100.0*var  # y is normalized to [0,1]. Sceled up by 100 to fit f.
	idx=int(scaly)  # find bin index corresponding to scaly
	f[idx] +=(idx+1-scaly)
	f[idx+1] +=(scaly-idx)

f=f/N  # normalization of f
plt.plot(f); plt.ylim(0,0.02); plt.show()
