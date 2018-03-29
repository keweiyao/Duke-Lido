#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import sys

mass = 1.3
T = 0.4
def norm_boltzmann(T):
	E = np.linspace(1, 10, 1000)*mass
	dfdE = np.exp(-E/T)*E*np.sqrt(E**2-mass**2)*4.*np.pi
	return np.sum(dfdE)*(E[1]-E[0])
def thermal_dfdE(T):
	E = np.linspace(1, 10., 1000)*mass
	dfdE = np.exp(-E/T)*E*np.sqrt(E**2-mass**2)*4.*np.pi
	return E, dfdE/norm_boltzmann(T)
		
ds = np.loadtxt(sys.argv[1]).reshape(50,10000,8)
x1 = ds[:,:,:4]
p1 = ds[:,:,4:]

ds = np.loadtxt(sys.argv[2]).reshape(50,10000,8)
x1 = ds[:,:,:4]
p2 = ds[:,:,4:]

E0 = p1[0,0,0]
x0, y0 = thermal_dfdE(T)
for it, (ip1, ip2) in enumerate(zip(p1, p2)):
	DE = E0-(ip1.T[0]).mean()
	DE2 = E0-(ip2.T[0]).mean()
	plt.plot(it*0.1, DE/E0, 'ro')
	plt.plot(it*0.1, DE2/E0, 'bD')
	#print(ip.shape)
	#plt.clf()
	#plt.hist(ip2.T[0], 100, histtype='step', normed=True)
	#plt.plot(x0, y0, 'r-')
	#plt.pause(.2)
plt.show()
	


