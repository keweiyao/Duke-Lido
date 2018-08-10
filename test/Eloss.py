#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import sys, h5py, os
from event import probe_run, Bjorken_run
import subprocess as sp

alphas_list = [0.075, 0.15, 0.3, 0.6]
T_list = [0.2, 0.4, 0.6]
E_list = np.exp(np.linspace(np.log(5), np.log(200), 20))
print(alphas_list)
print(T_list)
print(E_list)
results = np.zeros([len(alphas_list), len(T_list), len(E_list), 100])
dL_output = np.zeros([len(alphas_list), len(T_list), len(E_list)])
for i, alphas in enumerate(alphas_list):
	g = np.sqrt(alphas*4*np.pi)
	for j, T in enumerate(T_list):
		for k, E in enumerate(E_list):			
			approx_qhat = g**2*(1.5*g**2*T**2)*T*np.log(10/(1.5*g**2*T**2))/3./np.pi
			Lc = np.sqrt(E/approx_qhat)/5.026
			dL = Lc/500
			Lmax = Lc*2
			print("-----", Lmax, "-----")
			Np = int(1000/alphas)
			mode = 'old' if os.path.exists('table.h5') else 'new'
			if j==0 and k==0:
				mode = 'new'
			dE = probe_run(1.3, E, T, dt=dL, Nsteps=int(Lmax/dL), Nparticles=Np, mode=mode, mu=0.6, alphafix=-1)
			dL_output[i,j,k] = dL*10
			results[i,j,k] = dE[::10]
			
			
			
with h5py.File("Eloss.h5", 'a') as f:
	name = 'coll-rad-old'
	if name in f:
		del f[name]
	gp = f.create_group(name)
	gp.attrs.create('alphas', alphas_list)
	gp.attrs.create('T', T_list)
	gp.attrs.create('E', E_list)
	gp.create_dataset("dE", data=results)

	if 'dL' in f:
		del f['dL']
	f.create_dataset("dL", data=dL_output)


