import random
import numpy as np
import matplotlib.pyplot as plt

N=int(input("Number of steps="))
n=int(input("Number of drunken men="))

xvar=np.arange(N)
x=np.zeros(N)
rms=np.zeros(N)

for m in range(n):
	for i in range(N-1):
		x[i+1] = random.choice([1,-1])+x[i]  
	for i in range(N):
		rms[i]=rms[i]+x[i]**2 

rms = np.sqrt(rms/n)

plt.plot(xvar,rms,xvar,np.sqrt(xvar)); plt.show()
