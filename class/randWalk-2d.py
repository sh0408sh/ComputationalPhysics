import random
import numpy as np
import matplotlib.pyplot as plt

plt.ion()

N=int(input("Number of steps="))
x=np.zeros([N,2])
mn = -2.0*np.sqrt(1.0*N) 
mx =  2.0*np.sqrt(1.0*N) 

plt.figure(figsize=(10,10))

for i in range(N-1):
	r=random.choice([1,-1,2,-2])
	if r==1 :
		x[i+1,0]=x[i,0]+1
		x[i+1,1]=x[i,1]
	elif r==-1:	
		x[i+1,0]=x[i,0]-1
		x[i+1,1]=x[i,1]
	elif r==2:
		x[i+1,0]=x[i,0]
		x[i+1,1]=x[i,1]+1
	else: 
		x[i+1,0]=x[i,0]
		x[i+1,1]=x[i,1]-1
#	plt.plot(x[:i-1,0], x[:i-1,1],lw=1.0) 
#	plt.xlim(mn,mx)
#	plt.ylim(mn,mx)
#	plt.draw()
#	plt.pause(0.0001)
#	plt.clf()


plt.ioff()
plt.plot(x[:,0],x[:,1],'r',lw=1.0)
plt.xlim(mn,mx)
plt.ylim(mn,mx)
plt.show()
