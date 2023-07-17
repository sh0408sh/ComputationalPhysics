import numpy as np
import matplotlib.pyplot as plt

#Ntr=int(input("N-transient="))
Ntr=200
Nit=int(input("N-iteration="))
r1=float(input("r_start="))
r2=float(input("r_end="))

rs=np.arange(r1,r2,0.0001)

r=[]
x=[]
for var in rs:
	xn=0.5
	for var2 in range(Ntr):
		xn=var*xn*(1.0-xn) 

	for var2 in range(Nit):
		xn=var*xn*(1.0-xn) 
		r.append(var)
		x.append(xn)

plt.plot(r,x,'bo',ms=0.05); plt.show()

	
