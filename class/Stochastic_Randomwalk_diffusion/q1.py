import numpy as np
import matplotlib.pyplot as plt
import random

# Number of particles to use for simulation
N = int(1000)

# Size of the unit walk step.
dx = float(0.1)
dy = float(0.1)

# Max steps of simulation
smax = int(1000)

# Number of bins for density distribution plot
nBin = 10

# Calculation of theoretical parameters
D = 0.5 * dx**2  # suppose dt=1
tau = 1.0 / (D * np.pi * np.pi)

# Zero-initialized particle positions
x = np.zeros((N, 2))


# Construct the initial sine-distribution of 10000 particles
i = 0
while i < N:
    V = random.random()
    T = random.random()

    val = np.sin(np.pi * V)

    if T <= val:
        x[i,0] = V
        i += 1

i = 0
while i < N:
    V = random.random()
    T = random.random()
    val = np.sin(np.pi * V)
    if T <= val:
        x[i,1] = V
        i += 1

print(x)


# Create bin meshes for density distribution
f = np.zeros((nBin + 1, nBin + 1))
bnX = np.arange(0, nBin + 1)
bnX = bnX / (1.0 * nBin)  # x-value for the distr-plot
bnY = np.arange(0, nBin + 1)
bnY = bnY / (1.0 * nBin)  # y-value for the distr-plot
bnX, bnY = np.meshgrid(bnX,bnY)

forward = np.array([dx,0.])
backward = np.array([-dx,0.])
right = np.array([0.,dy])
left = np.array([0.,-dy])

fig=plt.figure()
#asx = fig.gca(projection='3d')
asx=fig.add_subplot(111,projection='3d')


for s in range(0):
    f = 0.0 * f  # empty the bin for new constructon of the distr. every step
    for i in range(N):
        if (0<=x[i,0] and x[i,0]<=1) and (0<=x[i,1] and x[i,1]<=1) :
            ax = x[i,0] * (1.0 * nBin)  # normalize the potision to bin
            ay = x[i,1] * (1.0 * nBin)
            mx = int(ax)  # find the bin index
            my = int(ay)
            ax -= mx  # find the relative displacement in bin
            ay -= my

            f[mx,my] += (1.0 - ax)*(1.0 -ay)
            f[mx, my+1] += (1.0 - ax)*ay
            f[mx+1, my] += ax*(1.0-ay)
            f[mx+1,my+1] +=  ax*ay

            x[i] = x[i] + random.choice([forward,backward,right, left])  # random walk
    f = (1.0 * nBin) * f  # normalization factor for distri.
    #z=0.5*np.pi*N*np.sin(np.pi*bnX)*np.sin(np.pi*bnY)*np.exp(-(1.0*s)/tau)
    #fl = lambda x,y : f[x , y]
    #asx.plot_surface(bnX, bnY, z)
    # asx.plot_surface(bnX, bnY, f)
    # asx.set_xlabel("x")
    # asx.set_ylabel("y")
    # asx.set_zlabel("f")
    # plt.draw(); plt.pause(0.01); plt.clf()
    asx.clear()
    # surf = ax.plot_surface(X, Y, z, cmap='gray')
    surf = asx.plot_surface(bnX, bnY, f, cmap='coolwarm', alpha=0.1)
    plt.pause(0.01)

plt.show()