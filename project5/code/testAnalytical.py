import numpy as np
import matplotlib.pyplot as plt

#%% Test analytical solution

dt = .1
dx = .01
T = .1
Nt = (int) (round(T/dt))
Nx = (int) (round(1./dx+1))
uSol = []
x = np.linspace(0,1,Nx)

dt = .01


for j in xrange(Nx):
    u = 0
    for k in xrange(1,Nx*100):
        u += np.exp(-(k*np.pi)**2*dt)*2./(k*np.pi)*(-1)**k*np.sin(k*np.pi*j*dx)
    u += j*dx
    uSol.append(u)

fig, ax = plt.subplots()
ax.plot(x, uSol)