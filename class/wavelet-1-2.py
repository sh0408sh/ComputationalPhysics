import pywt  # To use WT in python
import numpy as np
import matplotlib.pyplot as plt

# Get tmax and dt. 
tmax = float(input("tmax="))
dt = float(input("dt="))
N=tmax/dt   # num. samples (i.e. num. data points) 

# time array, going from 0 thru tmax, with interval dt
t = np.arange(0,N);  t *= dt    

# Create a chirped signal. Note: When you have a signal sin[phi(t)], how can 
# you determine its frequency in general? Here phi(t) is a time-dependent
# phase of the sine function. You may be familiar with phi(t)=w*t, where w
# is the angular frequency. When phi is not a linear function of t, what is
# the best way to determine its frequency? The answer is time derivative of 
# phi. For linear case, phi=w*t, d phi/dt=w, i.e. you get w from derivative.
# Extending this to the general case, you can define w := d phi /dt.
# In the signal below, phi=10*pi*t**2. Taking derivative, w=20*pi*t, i.e.
# frequency is not a constant but a function of time. In this way, you can
# create a time-varying-frequency.  
# Note that in this example, frequency goes from 0 thru 20*pi*tmax.
x = np.sin(10*np.pi * t**2 )
plt.plot(t,x) ; plt.show()

# Calculate continuous wavelet transform
# The 2nd argument of cwt() is an array of scale factor 'a' (see note 9-2). 
# When the sampled signal is transferred to cwt(), the dt information 
# is not transferred together, i.e. cwt() does not know actual sampling 
# period. Hence cwt() assumes dt=1. To convert the result to the actual 
# frequency, you have to divide the result by dt. Because the shortest period
# of the signal is 1, choosing 1 as the minimum value of 'a' is reasonable,
# though not always good. (remember, 'a' roughly corresponds to the width of 
# the wavelet). Hence it is usually good to make 'a'-array start from 1, 
# e.g. a=arange(1,50), but it is not always the optimum.  
#
# Another problem is that which frequency a given value of 'a' represents best
# is not readily known. Actually it depends on the type of the wavelet. 
# Try comparing the results of the following two lines (remove # temporarily).
#print(pywt.scale2frequency('cmor1.5-1.0',[1,2,3,4,5]))
#print(pywt.scale2frequency('mexh',[1,2,3,4,5]))
#quit()
# You may get the following results (tmax,dt do not matter to get this result)
#[1.         0.5        0.33333333 0.25       0.2       ]
#[0.25       0.125      0.08333333 0.0625     0.05      ]
# The function scale2frequency() returns an array of frequencies (with respect 
# to sampling period 1), that correspond to the values of 'a'-array (or list). 
# When 'cmor1.5-1.0' wavelet is chosen, a=1 correspond to frequency=1, while
# it is 0.25 for 'mexh' wavelet. Overall, the a-values 1,2,3,4,5 covers 4 times
# higher frequency band in 'cmor1.5-1.0' wavelet.  As mentioned, if you divide
# the result by dt, you get the actual frequencies.  
# Hence it is a useful method to pre-check which frequency band is examined
# by a given array of scale parameter 'a' for a given wavelet type.

# Determine the scale parameter array ('a'-array)
scl = np.arange(4,100,0.2)
#print(pywt.scale2frequency('cmor1.5-1.0',scl)/dt)  # check the frequency range

# Take the wavelet transform
# The 2nd return 'freqs' has the frequency information that is exactly the same
# as obtained by scale2frequency('...',scale)/dt. Here dt should be transferred
# to 'sampling_period' option of cwt() function. Note that 'sampling_period'
# does not affect the scale array. It is optional.
# In the following, you can compare WT by 'cmor1.5-1.0' and by 'mexh'.
# You can try different wavelets.
coef, freqs = pywt.cwt(x, scl, 'cmor1.5-1.0',sampling_period=dt)
#coef, freqs = pywt.cwt(x, scl, 'mexh',sampling_period=dt)

# Get the absolute square of the wt results. Squaring gives higher contrast
amp = np.abs(coef)**2

# pcolor() plots color-map for non-evenly distributed data points
plt.pcolor(t, freqs, amp); plt.show()

