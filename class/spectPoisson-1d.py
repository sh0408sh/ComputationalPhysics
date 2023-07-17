import numpy as np
from scipy.fftpack import fft, ifft
import matplotlib.pyplot as plt

xmax=float(input("xmax="))
N=int(input("Num. meshes="))  # N should be a power of 2 to use FFT 
Nhf=int(N/2)   # N/2 is to be used for arranging the k-array

dx=xmax/(1.0*N)
x=np.arange(0.0,N)
x *= dx  # get the x-domain in real scale

# arrange k-array
j= np.arange(0,N)
k= np.zeros(N)
k[:Nhf]=2.0*np.pi*j[:Nhf]/(1.0*N*dx)  # 1st half of k: [0, fc/2-df]
k[Nhf:]=2.0*np.pi*(j[Nhf:]-N)/(1.0*N*dx) # 2nd half of k: [-fc/2,-df]

# sine source
y=np.sin(2.0*np.pi*x/xmax)
plt.plot(x,y);  plt.show()  # visualize the charge source

# FFT the source
yh=fft(y)
yh[1:] /= (k[1:])**2  # Calculate yh for k!=0. eps0 is set to 1
yh[0]=0.0  # k=0 is the dc-component of the potential. It does not affect 
           # E-field as its derivative vanishes. Hence yh[0] is set to zero 
           # without losing information.  

# Finally get the potential from inverse FT of yh
z=ifft(yh)

plt.plot(x,z.real); plt.show()
