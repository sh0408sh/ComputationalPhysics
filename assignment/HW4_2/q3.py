import numpy as np
import matplotlib.pyplot as plt
import diffusion-1d-fcts-HH-Copy.1 as df

# Get the parameters
xmax=float(2) # Domain size
dx=float(0.1)   # Mesh size
dt=float(0.01)   # Time step
D=float(0.5)   # Diffusion coefficient
imax=int(500)  # max step to run

print(df.diffusion_eq(xmax,dx,dt,D,imax))