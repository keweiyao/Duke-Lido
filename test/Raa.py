#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import event, h5py
import fonll.fonll as fonll
from JetCalc.ExpCut import cuts

bins = np.concatenate([[cuts['CMS']['Raa']['pTbins'][0,0]],
					   cuts['CMS']['Raa']['pTbins'][:,1]])
mass = 1.3
sp = fonll.FONLL()

def norm_boltzmann(T):
	E = np.linspace(1, 20, 1000)*mass
	dfdE = np.exp(-E/T)*E*np.sqrt(E**2-mass**2)*4.*np.pi
	return np.sum(dfdE)*(E[1]-E[0])

def thermal_dfdE(T):
	E = np.linspace(1, 5., 1000)*mass
	dfdE = np.exp(-E/T)*E*np.sqrt(E**2-mass**2)*4.*np.pi
	return E, dfdE/norm_boltzmann(T)


e1 = event.event(medium={'type':'dynamic', 'hydrofile':"./JetData.h5"},
				 preeq={'type':'dynamic', 'hydrofile':"./FreeStream.h5"},
				 LGV={'A':0., 'B':2.}
				)
TAB = h5py.File('initial.hdf','r')['event_0/Ncoll_density'].value.T
MD = h5py.File('initial.hdf','r')['event_0/matter_density'].value.T
plt.plot(TAB[120])
plt.plot(MD[120])
print((TAB).sum()*0.135**2)
plt.show()
e1.initialize_HQ(30000, init_flags={'type':'A+B', 'pTmin':1, 'pTmax':150,
									'ymax':1.0, 'TAB':TAB, 'dxy':0.135, 'b':3.048})

p0, x0 = e1.HQ_hist()
pT0 = np.sqrt(p0.T[1]**2 + p0.T[2]**2)
w0 = sp.interp('EPPS', 'pp', '5020', 'c', pT0, np.zeros_like(pT0))
w = sp.interp('EPPS', 'PbPb', '5020', 'c', pT0, np.zeros_like(pT0))
H0, b0 = np.histogram(pT0, bins, weights=w0[:,0],range=[0,100])
H1, b0 = np.histogram(pT0, bins, weights=w[:,0],range=[0,100])
bm = (b0[:-1]+b0[1:])/2.
def plot():
	p, x = e1.HQ_hist()
	plt.clf()
	#plt.plot(E, dfdE, 'r-')
	pT = np.sqrt(p.T[1]**2 + p.T[2]**2)
	H, b = np.histogram(pT, bins, weights=w[:,0],range=[0,100])
	plt.plot(bm, H/H0, 'ro-')
	plt.plot(bm, H1/H0, 'r--')
	plt.semilogx()
	plt.ylim(0,1.2)
	plt.pause(0.2)
# loop through preeq:
for i in range(600):
	status = e1.perform_fs_step()
	print('{:1.2f}'.format(e1.sys_time()), '[fm/c]')
	if not status:
		break
	if i%1==0:
		plot()

for i in range(600):
	status = e1.perform_hydro_step()
	print('{:1.2f}'.format(e1.sys_time()), '[fm/c]')
	if not status:
		break
	if i%1==0:
		plot()

plt.show()

