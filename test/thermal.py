#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import event, h5py, os

T = 0.4
mass = 1.3

os.environ['XDG_DATA_HOME'] = './'
def norm_boltzmann(T):
	E = np.linspace(1, 20, 1000)*mass
	dfdE = np.exp(-E/T)*E*np.sqrt(E**2-mass**2)*4.*np.pi
	return np.sum(dfdE)*(E[1]-E[0])

def thermal_dfdE(T):
	E = np.linspace(1, 5., 1000)*mass
	dfdE = np.exp(-E/T)*E*np.sqrt(E**2-mass**2)*4.*np.pi
	return E, dfdE/norm_boltzmann(T)

x0, y0 = thermal_dfdE(T)

e1 = event.event(medium={'type':'static', 'static_dt':0.025},
				 LBT={'mu': 1.0}
				)

e1.initialize_HQ(N_charm=10000, N_bottom=1000,
				init_flags={'type':'probe', 'E0':10.})
tt = 0
with h5py.File("thermalization.hdf5", 'a') as f:
	name = 'coll-rad-abs'
	if name in f:
		del f[name]
	gp = f.create_group(name)
	for i in range(3000):
		status = e1.perform_hydro_step()
		print('{:1.2f}'.format(e1.sys_time()), '[fm/c]')
		if not status:
			break
		if i%40==0:
			p, x = e1.HQ_hist(4)
			gp.create_dataset("{:d}".format(tt), data=p.T)
			tt += 1


plt.show()
