import random
import scipy.signal as sg
import numpy as np
import matplotlib.pyplot as plt

# this is a signal with random noise
x=np.arange(100)
y0=np.sin(10.0*np.pi*x/100)
y=np.sin(10.0*np.pi*x/100)

for i in x:
	y[i] += 2.0*(random.random()-0.5)

plt.plot(x,y); plt.show()
##############################3


# Write your own sav-gol filter here and plot it along with the original signal


