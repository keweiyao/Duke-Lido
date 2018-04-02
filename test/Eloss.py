#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import sys, h5py, os
from event import probe_run, rate_run

dt = 0.1
Nsteps = 50
Np = 10000
L = np.arange(Nsteps)*dt
E = np.exp(np.linspace(np.log(2), np.log(100), 20))
f = h5py.File("Eloss.h5", 'a')
name = 'rad-abs'
if name in f:
	del f[name]
gp = f.create_group(name)
gp.attrs.create('T', [0.2, 0.4, 0.6])
gp.attrs.create('E', E)
gp.attrs.create('L', L)
for iT, T in enumerate([0.2, 0.4, 0.6]):
	for iE, E0 in enumerate(E):
		mode = 'old' if os.path.exists('table.h5') else 'new'
		R = probe_run(E0, T, dt=dt/3, Nsteps=Nsteps*3, Nparticles=Np, mode=mode)
		R = np.array(R)
		R = np.array(R[::3])
		gp.create_dataset("{}/{}".format(iT, iE), data=R)
f.close()
