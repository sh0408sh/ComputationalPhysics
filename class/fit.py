import numpy as np
import matplotlib.pyplot as plt

x=np.linspace(0,9,num=10)
y=np.array([1.1,2.2,3.1,2.9,4.1,6.7,9.8,11.1,13.2,15.3])

p=np.polyfit(x,y,1)
f=np.poly1d(p)

plt.plot(x,y,'o')
plt.plot(x,f(x))
plt.show()
