#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import sys, h5py, os
from event import probe_run

dt = 0.005
Nsteps = 1000
Np = 2000
L = np.arange(Nsteps)*dt
E = np.exp(np.linspace(np.log(2), np.log(100), 20))
f = h5py.File("Eloss-BDMPSZ.h5", 'a')
name = 'rad'
if name in f:
	del f[name]
gp = f.create_group(name)
gp.attrs.create('T', [0.2, 0.4, 0.6])
gp.attrs.create('E', E)
gp.attrs.create('L', L[::5])
for iT, T in enumerate([0.2, 0.4, 0.6]):
	for iE, E0 in enumerate(E):
		mode = 'old' if os.path.exists('table.h5') else 'new'
		dE = probe_run(E0, T, dt=dt, Nsteps=Nsteps, Nparticles=Np, mode=mode)
		gp.create_dataset("{}/{}".format(iT, iE), data=dE[::5])
f.close()
