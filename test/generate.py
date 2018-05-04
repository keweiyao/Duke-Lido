#!/usr/bin/env python3
import numpy as np
import subprocess
import h5py
import event
import os
from scipy.interpolate import interp2d
import matplotlib.pyplot as plt
os.environ['XDG_DATA_HOME'] = './'

def get_kappa(ptype, A, B):
	mass = 1.3 if ptype == 'c' else 4.2
	with h5py.File('./table.h5', 'r') as f:
		processes =['{}q2{}q'.format(ptype, ptype),
					'{}g2{}g'.format(ptype, ptype)]
		a = f['Boltzmann/'+processes[0]+'/rate/scalar'].attrs
		E = np.linspace(a['low-0'], a['high-0'], a['shape-0'])
		T = np.linspace(a['low-1'], a['high-1'], a['shape-1'])
		kappa_zT3 = np.zeros([len(E),len(T)])
		kappa_TT3 = np.zeros([len(E),len(T)])
		for process in processes:
			dNdt = f['Boltzmann/'+process+'/rate/scalar/0'].value
			dpzdt = f['Boltzmann/'+process+'/rate/vector/3'].value
			dpx2dt = f['Boltzmann/'+process+'/rate/tensor/5'].value
			dpz2dt = f['Boltzmann/'+process+'/rate/tensor/15'].value
			kappa_zT3 += (dpz2dt - dpzdt**2/dNdt)/T**3
			kappa_TT3 += dpx2dt/T**3
		return E, T, kappa_zT3, kappa_TT3, A+B/np.outer(E, T)
"""
mu = 0.893
A = 0.362
B = 0.764
event.event(medium={'type':'static', 'static_dt':0.05}, LBT={'mu':mu})
E, T, kz, kT, kd = get_kappa('c', A, B)
print(E.shape, T.shape, kT.shape)
qhat_T3 = interp2d(E, T, (kT.T+kd.T)*2.)
qhat = qhat_T3((10**2+1.3**2)**0.5, T)[:,0]
plt.plot(T/0.154, qhat )
plt.xlim(0.5, 4.0)
plt.xlabel(r"$T/T_c$",fontsize=15)
plt.ylabel(r"$2\pi T D_s$",fontsize=15)
plt.ylim(0,20)
plt.show()
with open("True-qhat.txt", 'w') as f:
	print("T/Tc\tqhat/T^3", file=f)
	for x, y in zip(T, qhat):
		print(x, '\t', y, file=f)


"""
with h5py.File('kappa.h5', 'a') as f:
	for nPDF in ['EPPS','nCTEQ']:
		mu, A, B = np.loadtxt(nPDF+'-sample-parameter.txt').T
		if nPDF in f:
			del f[nPDF]
		gpPDF = f.require_group(nPDF)
		for i, (m, a, b) in enumerate(zip(mu, A, B)):
			print(i, m, a, b)
			gp_name = "{:d}".format(i)
			gp = gpPDF.require_group(gp_name)
			gp.attrs.create('mu', mu)
			gp.attrs.create('A', a)
			gp.attrs.create('B', b)

			subprocess.call("rm table.h5", shell=True)
			event.event(medium={'type':'static', 'static_dt':0.05}, LBT={'mu':m})
			for ptype in ['c', 'b']:
				E, T, kz, kt, kd = get_kappa(ptype, a, b)
				gp.create_dataset('{}/kd'.format(ptype), data=kd)
				gp.create_dataset('{}/kt'.format(ptype), data=kt)
				gp.create_dataset('{}/kz'.format(ptype), data=kz)
			gp.attrs.create('E', E)
			gp.attrs.create('T', T)



