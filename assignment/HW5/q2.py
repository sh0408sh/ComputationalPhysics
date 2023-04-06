import numpy as np
import matplotlib.pyplot as plt
import random

# Number of particles to use for simulation
# N=int(input("Number of particles="))
N = 10000

# Size of the unit walk step.
# dx=float(input("Size of the walk step="))
dx = 1e-1

# Max steps of simulation
# smax=int(input("Max steps="))
smax = 1000

# Number of bins for density distribution plot
nBin = 10

# Calculation of theoretical parameters
D = 0.5 * (dx) ** 2  # suppose dt=1
tau = 1.0 / (D * np.pi * np.pi)
print( "theoretical characteristic time :", tau)

# x=np.zeros(N)  # Zero-initialized particle positions
# Construct the initial sine-distribution of 10000 particles
# i=0
# P = np.sin(np.pi*x)
# while i<N:
# 	V=random.random()
# 	T=random.random()
# 	val=np.sin(np.pi*V)
# 	if T<= val:
# 		x[i]=V
# 		i+=1

t = np.arange(N) / N
# P = np.sin(np.pi * t)
# P /= np.sum(P)
x = t*0
y = t*0
#
# # Create bin meshes for density distribution
f = np.zeros((nBin + 1, nBin + 1))
bnX = np.arange(0, nBin + 1)
bnX = bnX / (1.0 * nBin)  # x-value for the distr-plot

i=0
while i< N:
    V=random.random()
    T=random.random()
    val=np.sin(np.pi*V)
    if T<= val:
        x[i] = V
        i+=1
i=0
while i< N:
    V=random.random()
    T=random.random()
    val=np.sin(np.pi*V)
    if T<= val:
        y[i] = V
        i+=1

X, Y = np.meshgrid(bnX, bnX)
fig = plt.figure()
ax1=fig.add_subplot(121,projection='3d')
ax2=fig.add_subplot(122,projection='3d')

# Main loop
for s in range(smax):
    f = 0.0 * f  # empty the bin for new constructon of the distr. every step
    for i in range(N):
        if 0 <= x[i] and x[i] < 1 and 0 <= y[i] and y[i] < 1:
            nx = x[i] * (1.0 * nBin)  # normalize the potision to bin
            m = int(nx)  # find the bin index
            nx -= m  # find the relative displacement in bin
            ny = y[i] * (1.0 * nBin)
            n = int(ny)
            ny -= n

            f[m, n] += (1.0 - nx) * (1.0 - ny)
            f[m + 1, n] += nx * (1.0 - ny)
            f[m, n + 1] += (1.0 - nx) * ny
            f[m + 1, n + 1] += nx * ny

            r = random.choice([0,1])
            if r == 0 :
                x[i] = x[i] + random.choice([dx, -dx])
            else :
                y[i] = y[i] + random.choice([dx, -dx])

            # x[i] = x[i] + random.choice([dx, -dx])
            # y[i] = y[i] + random.choice([dx, -dx])  # random walk

    f = 1.0 * nBin * 1.0 * nBin * f  # normalization factor for distri.

    if s == 0 :
        fsum_0 = np.sum(f)
        print(fsum_0)

    if np.sum(f)/fsum_0 <= 1/np.e :
        tau = s
        print(np.sum(f))
        print(np.sum(f)/fsum_0,1/np.e )
        print( "measured characteristic time:",tau)
        break
    plt.ylim(0,1.2*np.pi*N/2.0)
    z = (np.pi / 2) ** 2 * N * np.sin(np.pi * X) * np.sin(np.pi * Y) * np.exp(-(1.0 * s) / tau)
    if s%20 == 0 :
        ax1.clear()
        ax2.clear()
        ax1.set_zlim(0, N * 2.5)
        ax2.set_zlim(0, N * 2.5)
        surf1 = ax1.plot_surface(X, Y, z, cmap='gray')
        surf2 = ax2.plot_surface(X, Y, f, cmap='coolwarm', alpha = 0.8)
        plt.pause(1)

# print(np.sum(f))
# print(np.sum(z))
plt.show()