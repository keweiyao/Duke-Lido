#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import event
import sys, os, h5py

os.environ['XDG_DATA_HOME'] = './'

# https://arxiv.org/pdf/1205.2396.pdf
initial = {'c': lambda x: x/(x**2+2.1**2)**3.9,
		   'b': lambda x: x/(x**2+7.5**2)**4.9 }

e1 = event.event(LBT={'mu':1.0}, medium={'type':'static', 'static_dt': 0.05})
e1.initialize_HQ(N_charm=400000, N_bottom=400000, init_flags={"type":"Box", "pmax": 15.})

pc0, xc0 = e1.HQ_hist(4)
pb0, xb0 = e1.HQ_hist(5)
for i in range(60):
	print(i, "steps")
	e1.perform_hydro_step({"Temp":0.3, "Vx":0.0, "Vy":0.0, "Vz":0.})

p1, x1 = e1.HQ_hist(4)
pTc0 = np.sqrt(np.sum(pc0.T[1:]**2, axis=0))
pTc1 = np.sqrt(np.sum(p1.T[1:]**2, axis=0))
wc = initial['c'](pTc0)

N=61
p1, x1 = e1.HQ_hist(5)
pTb0 = np.sqrt(np.sum(pb0.T[1:]**2, axis=0))
pTb1 = np.sqrt(np.sum(p1.T[1:]**2, axis=0))
wb = initial['b'](pTb0)
H0, x0 = np.histogram(pTc0, N, weights=wc, normed=True, range=[0,15])
H1, x0 = np.histogram(pTc1, N, weights=wc, normed=True, range=[0,15])
x0 = (x0[1:]+x0[:-1])/2.
raac = H1/H0

H0, x0 = np.histogram(pTb0, N, weights=wb, normed=True, range=[0,15])
H1, x0 = np.histogram(pTb1, N, weights=wb, normed=True, range=[0,15])
x0 = (x0[1:]+x0[:-1])/2.
raab = H1/H0

plt.plot(x0, raac, 'r-')
plt.plot(x0, raab, 'b-')
plt.ylim(0,1.5)
plt.show()
with open("BoxRaa-mu1-charm-coll-rad-abs-s.dat",'w') as f:
	for ix, iR in zip(x0, raac):
		print(ix, iR, file=f)
with open("BoxRaa-mu1-bottom-coll-rad-abs-s.dat",'w') as f:
	for ix, iR in zip(x0, raab):
		print(ix, iR, file=f)
