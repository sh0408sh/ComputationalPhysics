import random
import numpy as np

N=int(input("number of random numbers="))

sum=0
for i in range(N):
	r1=random.random()  # x-coord
	r2=random.random()  # y-coord
	r=np.sqrt(r1**2 + r2**2)
	if r<=1.0:  # accept when distance to (r1,r2) is less than the radius
		sum +=1

sum /= (1.0*N)
print(4.0*sum)  # print approximated pi value
