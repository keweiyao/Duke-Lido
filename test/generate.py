#!/usr/bin/env python3
import xml.etree.cElementTree as ET
import numpy as np
import subprocess as sp
import h5py

nPDF = 'EPS09'
parameters = np.loadtxt(nPDF+'-sample-parameter.txt').T
mu, A, B = parameters
fo = h5py.File('qhat_pQCD', 'a')
for i, m in enumerate(mu):
	tree = ET.parse('settings.xml')
	root = tree.getroot()
	root[0][0].text=str(m)
	tree.write('settings.xml')
	sp.call("hybrid", shell=True)
	f = h5py.File('./qhat.h5', 'r')
	gp = fo.create_group(nPDF+'/{}/'.format(i))
	
	for process,c in zip(['Qq2Qq', 'Qg2Qg'], ['r','g']):
		a = f[process+'/rate/scalar'].attrs
		E = np.linspace(a['low-0'], a['high-0'], a['shape-0'])
		T = np.linspace(a['low-1'], a['high-1'], a['shape-1'])
		R = f[process+'/rate/scalar/0'].value
		Rz = f[process+'/rate/vector/3'].value
		Rxx = f[process+'/rate/tensor/5'].value
		Rzz = f[process+'/rate/tensor/15'].value
		Rzz = Rzz - Rz**2/R
		gp.create_dataset(process+'/qxx', data=Rxx)
		gp.create_dataset(process+'/qzz', data=Rzz)
	f.close()
fo.close()
