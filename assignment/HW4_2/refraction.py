import numpy as np
import matplotlib.pyplot as plt
import emitEMwave as ant
import emitEMwave_update as antu
import math as mt


def refraction(xmax, ymax, dx, dy, dt, f, D, smax, phi, apk, n1, n2, dur):

    dsav = 50
    a = dt / dx;
    b = dt / dy
    w = 2.0 * np.pi * f
    s = 0

    cntr = 0.6 * ymax
    # upper=int(0.5*(ymax+D)/dy)
    # lower=int(0.5*(ymax-D)/dy)
    upper = int(D / dy)
    lower = int(0)

    x = np.arange(0, xmax + dx, dx)
    y = np.arange(0, ymax + dy, dy)
    Nh = int(xmax / dx / 2)

    n = n2 / n1  # relative refractive idnex

    X, Y = np.meshgrid(x, y)
    Ex = 0 * X;
    Ey = 0 * X;
    Ez = 0 * X
    Bx = 0 * X;
    By = 0 * X;
    Bz = 0 * X
    while s < smax:
        # Ey[:,0]= np.exp(-(y-cntr)**2/(0.2*ymax)**2)*np.sin(w*s*dt) # emission
        # Ey[lower:upper,0]= np.sin(w*s*dt) # hole
        ant.emitEMwave(s * dt, Ey[lower:upper, 0], (dx, dy), 'ovwrt', 'p', phi, 0, f, 0, 0, 6, (1.0,), dur, apk)
        ant.emitEMwave(s * dt, Ez[lower:upper, 0], (dx, dy), 'sum', 's', phi, 0, f, 0, 0, 6, (1.0,), dur, apk)
        Bx[:-1, :-1] += -b * (Ez[1:, :-1] - Ez[:-1, :-1])
        By[1:-1, :-1] += a * (Ez[1:-1, 1:] - Ez[1:-1, :-1])
        Bz[:-1, :-1] += -a * (Ey[:-1, 1:] - Ey[:-1, :-1]) + b * (Ex[1:, :-1] - Ex[:-1, :-1])

        Ex[1:-1, :Nh] += 1.0 * b * (Bz[1:-1, :Nh] - Bz[:-2, :Nh])
        Ex[1:-1, Nh:-1] += 1.0 * b * (Bz[1:-1, Nh:-1] - Bz[:-2, Nh:-1]) / n ** 2

        Ey[:-1, 1:Nh] += -1.0 * a * (Bz[:-1, 1:Nh] - Bz[:-1, :Nh - 1])
        Ey[:-1, Nh:-1] += -1.00 * a * (Bz[:-1, Nh:-1] - Bz[:-1, Nh - 1:-2]) / n ** 2

        Ez[1:-1, 1:Nh] += 1.0 * (a * (By[1:-1, 1:Nh] - By[1:-1, :Nh - 1]) - b * (Bx[1:-1, 1:Nh] - Bx[:-2, 1:Nh]))
        Ez[1:-1, Nh:-1] += 1.00 * (a * (By[1:-1, Nh:-1] - By[1:-1, Nh - 1:-2]) - b * (Bx[1:-1, Nh:-1] - Bx[:-2, Nh:-1])) / n ** 2
        s += 1
        '''
        if s % dsav == 0:
            cs = plt.imshow(Ey**2); plt.colorbar(cs)
            #plt.clim(0,1)
            plt.draw(); plt.pause(0.01); plt.clf()

        '''



    """"
    extn = (0, xmax, 0, ymax)
    plt.subplot(1, 2, 1);
    cs = plt.imshow(Ey, origin='lower', extent=extn);
    plt.colorbar(cs)
    plt.subplot(1, 2, 2);
    cs = plt.imshow(Ez, origin='lower', extent=extn);
    plt.colorbar(cs)
    plt.show()
    """

    """"
    if n1 > n2 :
        #print(np.max(abs(Ey[:,-1])))
        #print(np.max(abs(Ey[:,Nh+int(Nh/2):])))
        m1 = np.max(abs(Ey[:,Nh+3:]))
        m2 = np.max(abs(Ey[:,Nh:]))
        return abs(m1-m2)
    else :
        print(np.sum(Ey[:,Nh]**2))
    """

    return np.sum(Ey[:, :Nh] ** 2)