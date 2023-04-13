import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft, ifft  # to use FFT and iFFT

# input parameters - tmax and dt
tmax=float(input("tmax="))
dt=float(input("dt="))

# Calculate number of data points from tmax and dt
N=int(tmax/dt)

# t-domain (corresponds to x for h(x) type) 
t= np.arange(0,N)  # t goes from 0 thru N-1
t = dt*t           # t goes from 0 thru tmax-dt

# f-domainrange which is to be used for x-axis of power spectrum plot
# f corresponds to k for fft of h(x) type data
df=1.0/(dt*N)      # frequency interval
f=np.arange(-N/2,N/2)  
f *=df

# sine source and 2nd harmonic
y=np.sin(16.0*np.pi*t)+ 0.5*np.sin(32*np.pi*t)  # frequency is 8 and 16
plt.plot(t,y);  plt.show()

# FT
yh=fft(y)
yh1=yh[0:int(N/2)]  # 1st half of fft'ed data; contains [0,f_c)
yh2=yh[int(N/2):]   # 2nd half of fft'ed data; contains [-f_c,0)
yhs=np.concatenate((yh2,yh1))  # put the 2nd half [-f_c,0) to the left of 1st half [0,f_c), so that yhs has frequency spectrum from [-f_c,f_c). 
pw = abs(yhs)**2   # take the absolute squre of the result to get 'power'
plt.plot(f,pw); plt.show()

# inverse FT
z=ifft(yh)  # reconstruct of the original signal by iFFT.

# the ifft yields complex values although complex part is 0 for real signal, but still complex type. plot() gets only 'real' type, i.e. float type. Hence z.real should be used. 
plt.plot(t,z.real); plt.show()  
