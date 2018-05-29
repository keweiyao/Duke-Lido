#!/usr/bin/env python3
import xml.etree.cElementTree as ET
import numpy as np
import matplotlib.pyplot as plt
import h5py

for Q in ['c', 'b']:
	f = h5py.File('./table-{}.h5'.format(Q), 'r')
	process = 'Qg2Qg'
	a = f[process+'/rate/scalar'].attrs
	E = np.linspace(a['low-0'], a['high-0'], a['shape-0'])
	T = np.linspace(a['low-1'], a['high-1'], a['shape-1'])
	R = f[process+'/rate/scalar/0'].value
	Rz = f[process+'/rate/vector/3'].value
	Rxx = f[process+'/rate/tensor/5'].value
	Rzz = f[process+'/rate/tensor/15'].value
	Rzz = Rzz - Rz**2/R
	D = 4*np.pi/(Rxx/T**3)
	plt.plot(E, D, color=Q)
plt.show()
