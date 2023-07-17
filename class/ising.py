import numpy as np
import matplotlib.pyplot as plt
import random

# Get some parameters
Nr=30
Nc=30
imax=100000
h=0.0
T=2.5

S=np.arange(0,Nr*Nc)
S = S.reshape(Nr,Nc)

for i in range(Nr):
	for j in range(Nc):
#		S[i,j] = random.choice([1,-1])  # hot start
		S[i,j] = 1.0   # cold start

E=M=0
for i in range(Nr):
	for j in range(Nc):
		if i==0:     sb = S[Nr-1,j] 
		else:        sb = S[i-1,j]

		if i==Nr-1:  st = S[0,j] 
		else:        st = S[i+1,j]
 
		if j==0:     sl = S[i,Nc-1]
		else:        sl = S[i,j-1]

		if j==Nc-1:  sr = S[i,0]
		else:        sr = S[i,j+1]

		E += -S[i,j]*(sb+st+sl+sr) - S[i,j]*h
		M += S[i,j]

# plot of Charge distribution 
cs=plt.imshow(S); plt.colorbar(cs); plt.show()

nSmpl=1
Eacc = E
Macc = M
Elist=[]
Mlist=[]
smplList=[]
while nSmpl< imax:

	i = random.choice(range(Nr))
	j = random.choice(range(Nc))

	if i==0:     sb = S[Nr-1,j] 
	else:        sb = S[i-1,j]

	if i==Nr-1:  st = S[0,j] 
	else:        st = S[i+1,j]

	if j==0:     sl = S[i,Nc-1]
	else:        sl = S[i,j-1]

	if j==Nc-1:  sr = S[i,0]
	else:        sr = S[i,j+1]

	dltM = -2.0*S[i,j]
	dltE = 2.0*S[i,j]*(sb+st+sl+sr)+ 2*h*S[i,j]
	prb=np.exp(-dltE/T)
	if prb >1:  prb=1

	if prb >= random.random():  
		S[i,j] *= -1
		E += dltE
		M += dltM
		Eacc += E
		Macc += M
		nSmpl +=1
		Eave = Eacc/nSmpl
		Mave = Macc/nSmpl
		Elist.append(Eave)
		Mlist.append(np.abs(Mave)/(Nr*Nc))
		#Mlist.append(M)
		smplList.append(nSmpl)
	
#	if l%100==0: 
#		cs=plt.imshow(S); plt.colorbar(cs); plt.draw(); plt.pause(0.1); plt.clf()

	if nSmpl%2000==0: 
		plt.plot(smplList,Mlist); 
		plt.ylim(0,1)
		plt.xlim(0,imax)
		plt.draw(); plt.pause(0.1); plt.clf()

cs=plt.imshow(S); plt.colorbar(cs); plt.show()


