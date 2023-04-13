# Create N random num obeying sin-function distribution
import random
import numpy as np
import matplotlib.pyplot as plt

N=int(input("Num points="))   # number of data points
nbin = int(input("Num bins in distribution="))  # get num bins used for distr

v=np.zeros(N)   # initialize an array of N elements by zero

i=0
while i<N:
	r1=random.random()  # gen a random num corr to x of sin func
	val=np.sin(np.pi*r1) # sin value at r1
	r2=random.random()  # gen another random, to compare with func value 
	if r2<=val:  # accept when r2 is below the functional value
		v[i]=r1
		i+=1  # Loop until num accepted points reaches N

plt.plot(v,'ro',ms=1); plt.show()

# make a plot of density distribution
f=np.zeros(nbin+1)
for i in range(N):
	scalv= v[i]*nbin
	m=int(scalv)  # find the index of the bin the rand num belongs to
	dx=scalv-m
	f[m+1] += dx;  f[m] +=(1.0-dx) 

mx = np.max(f); f = f/mx  # Normalize f by its max, so that f peak becomes 1
x=np.arange(0.0,1.0*nbin+1.0);  x *= (1.0/nbin)

z = np.sin(np.pi*x)
plt.plot(x,f,x,z); plt.ylim(0,1.2); plt.show()
