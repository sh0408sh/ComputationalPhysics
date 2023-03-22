import numpy as np
import matplotlib.pyplot as plt

xmax=float(4)
dx=float(0.01)
dt=float(0.005)
f=float(4)
smax=int(800)
n1 = float(1)
n2 = float(1.3)
dsav = int(50)

a=dt/dx
w=2.0*np.pi*f

x=np.arange(0,xmax+dx,dx)
c=int(0.5*xmax/dx)

Ey=0*x; Ez=0*x
By=0*x; Bz=0*x

Ei = 0
tmp = 0
s = 0

while s < smax:
	By[:-1] += a * (Ez[1:] - Ez[:-1])
	Bz[:-1] += -a * (Ey[1:] - Ey[:-1])

	Ey[1:c] += -a * (Bz[1:c] - Bz[0:c - 1]) / (n1 ** 2)
	Ey[c:-1] += -a * (Bz[c:-1] - Bz[c - 1:-2]) / (n2 ** 2)
	if w * s * dt < 4.05/2*np.pi:
		Ey[0] += dt * (np.sin(w * s * dt))/ (n1 ** 2)
		#Ey[0] =  (np.sin(w * s * dt))
		tmp = Ey[0]
		if Ei < tmp :
			Ei = tmp
	Ez[1:c] += a * (By[1:c] - By[0:c - 1]) / (n1 ** 2)
	Ez[c:-1] += a * (By[c:-1] - By[c - 1:-2]) / (n2 ** 2)

	#if max < Ey[c+3] : max = Ey[c+3]
	#if min > Ey[c-3] : min = Ey[c-3]

	if s % dsav == 0:
		plt.ylim(-0.12,0.12)  # set the ylimit of sub-panels
		plt.yticks(np.arange(-0.12,0.12, 0.02))  # yticks
		plt.plot(x, Ey)
		plt.draw()
		plt.pause(0.01)
		plt.clf()

	s += 1


#plt.plot(x,Ey)
plt.show()

Et = np.max(Ey[c:])
Er = abs(np.min(Ey[:c]))
print(f"Ei is {Ei}\nEt is {Et}\nEr is {Er}")

m_T = (n2/n1)*((Et/Ei)**2)
m_R = (Er/Ei)**2
print(f"measured T is {m_T}\nmeasured R is {m_R}\n T+R={m_T+m_R}")

c_T = 4*(n1*n2)/((n1+n2)**2)
c_R = ((n1-n2)/(n1+n2))**2
print(f"calculated T is {c_T}\ncalculated R is {c_R}\n T+R={c_T+c_R}")

err_T = (abs(m_T-c_T)/c_T)*100
err_R = (abs(m_R-c_R)/c_R)*100
print(f"Error of T is {err_T}%\nError of R is {err_R}%")