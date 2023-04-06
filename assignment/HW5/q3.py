import numpy as np
import matplotlib.pyplot as plt
import random

# Number of particles to use for simulation
N = int(10000)

# Size of the unit walk step.
dx = float(0.1)

# Max steps of simulation
smax = int(1000)

# Number of bins for density distribution plot
nBin = 10

# Calculation of theoretical parameters
D = 0.5 * (dx) ** 2  # suppose dt=1
tau = 1.0 / (D * np.pi * np.pi)

x = np.zeros(N)  # Zero-initialized particle positions

# Construct the initial sine-distribution of 10000 particles
i = 0
while i < N:
    V = random.random()
    T = random.random()
    val = np.sin(np.pi * V)
    if T <= val:
        x[i] = V
        i += 1

w = [x]*4

print(x)
# Create bin meshes for density distribution
f = np.zeros(nBin + 1)
fw = [f]*4
bnX = np.arange(0, nBin + 1)
bnX = bnX / (1.0 * nBin)  # x-value for the distr-plot

#standard deviation
sigma = float(0.2)
# Main loop
for s in range(smax):
    f = 0.0 * f  # empty the bin for new constructon of the distr. every step
    for i in range(N):
        if 0 <= x[i] and x[i] <= 1:
            y = x[i] * (1.0 * nBin)  # normalize the potision to bin
            m = int(y)  # find the bin index
            y -= m  # find the relative displacement in bin
            f[m] += (1.0 - y);
            f[m + 1] += y
            #x[i] = x[i] + random.choice([dx, -dx])
            x[i] = x[i] + np.random.normal(0,sigma,1)[0]  # random walk
    f = (1.0 * nBin) * f  # normalization factor for distri.

    if s == 0 :
        fsum_0 = np.sum(f)
        print(fsum_0)

    if np.sum(f)/fsum_0 <= 1/np.e :
        ct = s
        print(np.sum(f))
        print(np.sum(f)/fsum_0,1/np.e )
        print(ct)
        break

    plt.ylim(0, 1.2 * np.pi * N / 2.0)
    z = 0.5 * np.pi * N * np.sin(np.pi * bnX) * np.exp(-(1.0 * s) / tau)
    plt.plot(bnX, f, c='r')
    plt.plot(bnX, z, c='b')
    plt.draw()
    plt.pause(1)
    plt.clf()

plt.ylim(0, 1.2 * np.pi * N / 2.0)
plt.plot(bnX, f);
plt.show();
