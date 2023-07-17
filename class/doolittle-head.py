import numpy as np
import scipy.linalg as la
A=np.array([[2,4,5,-1,-3],
           [1,3,2,6,-2],
           [4,2,1,0,-3],
           [-3,0,6,1,2],
           [-1,1,2,4,2]])


n=5
u=np.zeros(n*n);  u=np.reshape(u,[n,n])
l=np.zeros(n*n);  l=np.reshape(l,[n,n])
B=np.zeros(n*n);  B=np.reshape(B,[n,n])




# Write your own doolittle decompositor here

########################






# From here, it is just to check LxU gives the original matrix
# matrix product of L and U
for i in range(n):
	for j in range(n):
		s=0; 
		for k in range(n): s += l[i,k]*u[k,j]
		B[i,j] = s

print(A)  # original matrix
print(l)  # lower triangular
print(u)  # upper triangular
print(B)  # Check if LU=A

