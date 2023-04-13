import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # for 3D plot

	
a=65539.0
b=0
c=2.0**31

N=10000
t=np.arange(0,N)
x=0.0*t
y=0.0*t
z=0.0*t


new_z=1
for var in t:
	new_x=(a*new_z)%c
	new_y=(a*new_x)%c
	new_z=(a*new_y)%c
	x[var]=new_x/c; y[var]=new_y/c; z[var]=new_z/c

fig=plt.figure()
ax=fig.add_subplot(111,projection='3d')

ax.scatter(x,y,z,s=0.1,c='r',marker='o',linewidths=None)

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

plt.show()

plt.plot(-9*x+6*y,z,'ro',ms=0.5); plt.show()
