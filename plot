#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import h5py
from numpy import *


f = h5py.File('./table.h5', 'r')
M = 1.3
"""
process = 'Boltzmann/cq2cqg'
a = f[process+'/xsection/scalar'].attrs
ss = np.linspace(a['low-0'], a['high-0'], a['shape-0'])
T = np.linspace(a['low-1'], a['high-1'], a['shape-1'])
dt = np.linspace(a['low-2'], a['high-2'], a['shape-2'])
print(dt)
R = f[process+'/xsection/scalar/0'].value
fmax = f[process+'/xsection/fmax/0'].value

for k, idt in enumerate(dt):
	plt.subplot(6, 5, k+1)
	for j, iT in enumerate(T):
		
		approx = (idt*iT)**2/((idt*iT)**2+2.)/iT**4
		plt.plot(ss, fmax[:,j,k]/approx,'r-')
plt.show()

for k, idt in enumerate(dt):
	plt.subplot(6, 5, k+1)
	for j, iT in enumerate(T):
		approx = (idt*iT)**2/((idt*iT)**2+1)/iT**3
		plt.plot(ss, R[:,j,k]/approx,'r-')

plt.show()



a = f[process+'/rate/scalar'].attrs
E = np.linspace(a['low-0'], a['high-0'], a['shape-0'])
T = np.linspace(a['low-1'], a['high-1'], a['shape-1'])
dt = np.linspace(a['low-2'], a['high-2'], a['shape-2'])

R = f[process+'/rate/scalar/0'].value
fmax = f[process+'/rate/fmax/0'].value

for k, idt in enumerate(dt):
	plt.subplot(6, 5, k+1)
	for j, iT in enumerate(T):
		
		approx = (idt*iT)**2/((idt*iT)**2+3)*iT
		plt.plot(E, fmax[:,j,k]/approx,'r-')
plt.show()

for k, idt in enumerate(dt):
	plt.subplot(6, 5, k+1)
	for j, iT in enumerate(T):
		approx = (idt*iT)**2/((idt*iT)**2+5)*iT
		plt.plot(E, R[:,j,k]/approx,'r-')

plt.show()

"""
process = 'Boltzmann/cqg2cq'

a = f[process+'/xsection/scalar'].attrs
sqrts = np.linspace(a['low-0'], a['high-0'], a['shape-0'])
T = np.linspace(a['low-1'], a['high-1'], a['shape-1'])
x = np.linspace(a['low-2'], a['high-2'], a['shape-2'])
y = np.linspace(a['low-3'], a['high-3'], a['shape-3'])
dt = np.linspace(a['low-4'], a['high-4'], a['shape-4'])
X = f[process+'/xsection/scalar/0'].value
fmax = f[process+'/xsection/fmax/0'].value

def xapprox(sqrts, T, x, y, dt):
	M = 1.3
	s = sqrts**2
	M2 = M*M
	s12 = x*(s-M2) + M2
	s1k = (y*(1-s12/s)*(1-M2/s12)+M2/s12)*s 
	sqrts12 = s12**0.5

	E1 = (s12+M2)/2./sqrts12
	p1 = (s12-M2)/2./sqrts12
	E3 = (s+M2)/2./sqrts
	p3 = (s-M2)/2./sqrts
	k = (s-s12)/2./sqrts12
	costhetak = (M2 + 2.*E1*k - s1k)/2./p1/k
	sinthetak = np.sqrt(1. - costhetak*costhetak)
	kt = k*sinthetak
	kt2 = kt*kt
	kz = k*costhetak
	x = (k+kz)/sqrts12
	xbar = (k+np.abs(kz))/sqrts12
	tauf = 2.*(1.-xbar)*k/(kt2+x*x*M2)
	#return 
	return (dt*T)**2/(1.+(dt*T)**2) / T**2 /sqrts, (dt*T)**2/(.2+(dt*T)**2) / T**2

for k, idt in enumerate(dt):
	plt.subplot(5, 4, k+1)
	for j, iT in enumerate(T):
		approx, amax = xapprox(sqrts, iT,x[1],y[5],idt)
		plt.plot(sqrts, fmax[:,j,1,5,k]/amax,'r-')
plt.show()



a = f[process+'/rate/scalar'].attrs
E = np.linspace(a['low-0'], a['high-0'], a['shape-0'])
T = np.linspace(a['low-1'], a['high-1'], a['shape-1'])
dt = np.linspace(a['low-2'], a['high-2'], a['shape-2'])
print(dt)
R = f[process+'/rate/scalar/0'].value
fmax = f[process+'/rate/fmax/0'].value

for k, idt in enumerate(dt):
	plt.subplot(5, 4, k+1)
	for j, iT in enumerate(T):
		approx = .00001/E*(idt*iT)**2/((idt*iT*.5)**2+1)
		plt.plot(E, fmax[:,j,k]/approx,'r-')
plt.show()

for k, idt in enumerate(dt):
	plt.subplot(5, 4, k+1)
	for j, iT in enumerate(T):
		approx = 1./E*idt**1*iT**3
		plt.plot(E, R[:,j,k]/approx,'r-')

plt.show()
			

