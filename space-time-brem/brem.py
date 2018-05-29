#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import math, os
import scipy.integrate as integrate
from scipy.optimize import fmin as minimize
from fonll.fonll import FONLL
# sampling a 1 D function with rejection method
def sample1d(f, low, high, args=(), fmax=1.0):
	x = np.random.rand()*(high-low) + low
	ncount = 0
	while f(x, *args)/fmax < np.random.rand() and ncount < 1000:
		x = np.random.rand()*(high-low) + low
		ncount += 1
	val = f(x, *args)/fmax
	if val > 1.0:
		print("warning f/fmax = {} > 1.0".format(val))
	if ncount >= 1000:
		print("warning too many tries")
	return x;

# fixed alpha_s
Lambda2 = 0.2**2
alpha0 = 4*np.pi/9.
alpha_max = alpha0/np.log(2.)
def alpha_s(Q2):
	if Q2 <= Lambda2:
		return alpha_max
	else:
		return alpha0/np.log(1.+Q2/Lambda2)

# dR/dx kernel
# u = t*Q2/2/E, v = t*M2/2/E, w = t*E2/2/E 
def kernel(x, u, v, w):
	if x<1e-3 or x >= 1.0-1e-3:
		return 0.0
	else:
		u = np.min([x**2*w, np.ones_like(x)*u], axis=0)
		return 4./3.*(1+(1-x)**2)/x*(np.exp(-(v*x)**2) - np.exp(-(u+v*x**2)**2/x**2))


def tab_rate(N=41):
	if os.path.exists("rate.dat"):
		rate_table = np.fromfile("rate.dat").reshape(N,N,N)
	else:
		rate_table = np.zeros([N,N,N])
		u = np.linspace(1e-2, 1e1, N)
		v = np.linspace(1e-2, 1e1, N)
		w = np.linspace(1e-2, 1e1, N)
		for i, eu in enumerate(u):
			for j, ev in enumerate(v):
				print(j)
				for k, ew in enumerate(w):
					res,err = integrate.quad(kernel, 0.,1.,args=(eu,ev,ew))
					rate_table[i,j,k] = res
		with open("tab.dat",'w') as f:
			rate_table.tofile(f)

	if os.path.exists("diff_rate_max.dat"):
		max_table = np.fromfile("diff_rate_max.dat").reshape(N,N,N)
	else:
		max_table = np.zeros([N,N,N])
		u = np.linspace(1e-2, 1e1, N)
		v = np.linspace(1e-2, 1e1, N)
		w = np.linspace(1e-2, 1e1, N)
		for i, eu in enumerate(u):
			for j, ev in enumerate(v):
				print(j)
				for k, ew in enumerate(w):
					#print("\t", k)
					x = minimize(lambda x: -kernel(x,eu,ev,ew), 
								0.5, xtol=0.001, disp=0)[0]
					max_table[i,j,k] = 1.1*kernel(x, eu, ev, ew)
		with open("diff_rate_max.dat",'w') as f:
			max_table.tofile(f)

	return rate_table, max_table, (10.0-0.01)/(N-1.), N

Rtab, Rmaxtab, d, N = tab_rate()

@np.vectorize
def rate(E, M, Q2, t):
	u = t*Q2/2./E
	v = t*M**2/2./E
	w = t*E/2.
	ru, iu = np.modf((u-0.01)/d)
	rv, iv = np.modf((v-0.01)/d)
	rw, iw = np.modf((w-0.01)/d)
	iu = np.max([np.min([iu,N-2]),0])
	iv = np.max([np.min([iv,N-2]),0])
	iw = np.max([np.min([iw,N-2]),0])
	ru = u - d*iu
	rv = v - d*iv
	rw = w - d*iw
	res = 0.
	rmax = 0.
	for r1, i1 in zip([1-ru/d,ru/d],[int(iu), int(iu)+1]):
		for r2, i2 in zip([1-rv/d,rv/d],[int(iv), int(iv)+1]):
			for r3, i3 in zip([1-rw/d,rw/d],[int(iw), int(iw)+1]):
				res += Rtab[i1,i2,i3]*r1*r2*r3
				rmax += Rmaxtab[i1,i2,i3]*r1*r2*r3
	return res*alpha_max/2./np.pi/t, rmax
# rate with alpha_max, rejection on alpha performed later
rate = np.vectorize(rate)

def IsRad(E, M, Q2, t, dt):
	u = t*Q2/2./E
	v = t*M**2/2./E
	w = t*E/2.
	r, rmax = rate(E, M, Q2, t)
	if r*dt > np.random.rand():
		x = sample1d(kernel, 0, 1, args=(u, v, w), fmax=1.5*rmax)
		#sample y = exp(-t/tau_f) ~ uniform
		ymax = np.exp(-(v*x)**2)
		u = np.min([x**2*w, np.ones_like(x)*u], axis=0)
		ymin = np.exp(-(u/x+v*x)**2)
		y = np.random.uniform(ymin, ymax)
		kt2 = np.sqrt(-np.log(y))*2*x*E/t - (x*M)**2
		# reject on alpha_s
		if alpha_s(kt2)/alpha_max < np.random.rand():
			return False, []
		else:
			k = x*E
			kz = np.sqrt(k**2 - kt2)
			kt = np.sqrt(kt2)
			phi = np.random.uniform(0, np.pi*2.)
			return True, [k, kt*np.cos(phi), kt*np.sin(phi), kz]
	else:
		return False, []

def simulate1(E0, M):
	dt = 0.1
	t = 0.01
	Q2 = (E0**2-M**2)
	X = []
	while Q2 > Lambda2 and t<100:
		t += dt
		status, kmu = IsRad(E0, M, Q2, t, dt)
		if status:
			kt2 = kmu[1]**2+kmu[2]**2
			X.append(kmu)
			p0 = np.sqrt(E0**2-M**2)	
			Q2 = kt2
			E0 = E0 - kmu[0]
	X = np.array(X)
	if X.size == 0:
		return [], [], [E0]
	kt = np.sqrt(X.T[1]**2 + X.T[2]**2)
	k = X.T[0]
	return kt, k, [E0]
KT, K, E = [], [], []

fonll = FONLL()
Np = 10000
M = 1.3
pT = 1+np.random.rand(Np)*150
E0 = np.sqrt(pT**2+M**2)
wFO = fonll.interp('FO-EPPS', 'pp', '2760', 'c', pT, np.zeros_like(pT))
wFONLL = fonll.interp('EPPS', 'pp', '2760', 'c', pT, np.zeros_like(pT))
for i in range(Np):
	if i%100 == 0:
		print(i)
	kt,k,e = simulate1(E0[i], M)
	KT = np.concatenate([KT, kt])
	K = np.concatenate([K, k])
	E = np.concatenate([E, e])
#plt.hist(K, 20, normed=True)
#plt.show()
plt.hist(KT, 20, normed=True)
plt.show()
NLO, bins, _ = plt.hist(E0, 25, weights=wFO, range=[0,150], normed=True, histtype='step',color='r')
FONLL, bins, _ = plt.hist(E0, 25, weights=wFONLL, range=[0,150], normed=True, histtype='step',color='k')
MC, bins, _ = plt.hist(E, 25, weights=wFO, range=[0,150], normed=True, histtype='step', color='b')
plt.semilogy()
plt.show()

b = (bins[1:]+bins[:-1])/2.
plt.plot(b, NLO/FONLL,'r-')
plt.plot(b, MC/FONLL,'b-')
plt.show()

