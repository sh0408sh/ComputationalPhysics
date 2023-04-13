import random
import numpy as np
import matplotlib.pyplot as plt

Trial=10000
N=int(input("Number of steps="))
x=np.zeros(N)

cnt=0
for t in range(Trial):
	for i in range(N-1):
		r=random.choice([1,-1])
		if r==1 :
			x[i+1]=x[i]+1
		else:	
			x[i+1]=x[i]-1
		if x[i+1]==0:
			print(t,i)
			cnt +=1
			break

print(cnt/Trial)		

