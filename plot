#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import h5py

f = h5py.File('./table.h5', 'r')
a = f['Qq2Qq/xsection/scalar'].attrs
ss = np.linspace(a['low-0'], a['high-0'], a['shape-0'])
T = np.linspace(a['low-1'], a['high-1'], a['shape-1'])

X = f['Qq2Qq/xsection/scalar/0'].value
dpzX = f['Qq2Qq/xsection/vector/3'].value
dpx2X = f['Qq2Qq/xsection/tensor/5'].value
dpy2X = f['Qq2Qq/xsection/tensor/10'].value
dpz2X = f['Qq2Qq/xsection/tensor/15'].value

plt.plot(ss, X[:,5],'rD')

#plt.plot(ss, dpx2X[:,5]/X[:,5],'r-')
#plt.plot(ss, dpz2X[:,5]/X[:,5]-dpzX[:,5]**2/X[:,5]**2,'r--')

#plt.plot(T, dpx2X[10]/X[10]/T/T,'r-')
#plt.plot(T, (dpz2X[15]/X[15]-dpzX[15]**2/X[15]**2)/T/T,'r--')

a = f['Qg2Qg/xsection/scalar'].attrs
ss = np.linspace(a['low-0'], a['high-0'], a['shape-0'])
T = np.linspace(a['low-1'], a['high-1'], a['shape-1'])

X = f['Qg2Qg/xsection/scalar/0'].value
dpzX = f['Qg2Qg/xsection/vector/3'].value
dpx2X = f['Qg2Qg/xsection/tensor/5'].value
dpy2X = f['Qg2Qg/xsection/tensor/10'].value
dpz2X = f['Qg2Qg/xsection/tensor/15'].value

plt.plot(ss, X[:,5],'bD')

#plt.plot(ss, dpx2X[:,5]/X[:,5],'b-')
#plt.plot(ss, dpz2X[:,5]/X[:,5]-dpzX[:,5]**2/X[:,5]**2,'b--')

#plt.plot(T, dpx2X[10]/X[10]/T/T,'b-')
#plt.plot(T, (dpz2X[15]/X[15]-dpzX[15]**2/X[15]**2)/T/T,'b--')
plt.show()

