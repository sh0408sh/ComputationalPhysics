import numpy as np
import matplotlib.pyplot as plt

xmax = float(4)
ymax = float(4)
dx = float(0.01)
dy = float(0.01)
dt = float(0.005)
f = float(2)
D = int(50)
smax = int(500)
dsav = int(100)

# Simulation parameters
alx=dt/dx   # alpha parameter
aly=dt/dy
ax=alx**2
ay=aly**2
b=2*(1-ax-ay)

w=2.0*np.pi*f

# Mesh
x=np.arange(0,xmax+dx,dx)
y=np.arange(0,ymax+dy,dy)
X,Y=np.meshgrid(x,y)

cntry=int(0.5*ymax/dy)
cntrx=int(0.5*xmax/dx)

u1= 0*X; ubuf1= 0*X; uo1= 0*X

s = 2
while s <= smax:
    u1[cntrx, cntry-D] = np.sin(w * s * dt); u1[cntrx, cntry+D] = np.sin(w * s * dt)
    ubuf1[1:-1, 1:-1] = u1[1:-1, 1:-1]
    u1[1:-1, 1:-1] = b * u1[1:-1, 1:-1] + ax * (u1[2:, 1:-1] + u1[0:-2, 1:-1]) + ay * (
                u1[1:-1, 2:] + u1[1:-1, 0:-2]) - uo1[1:-1, 1:-1]
    uo1[1:-1, 1:-1] = ubuf1[1:-1, 1:-1]

    if s % dsav == 0: cs = plt.imshow(u1); plt.colorbar(cs); plt.clim(-1, 1); plt.draw(); plt.pause(
        0.01); plt.clf()
    s += 1

cs=plt.pcolormesh(X,Y,u1); plt.colorbar(cs); plt.show()