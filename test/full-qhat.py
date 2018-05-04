#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import event, h5py, os

os.environ['XDG_DATA_HOME'] = './'

def decomp(p0, p): # put in three vector
	p0abs = np.sqrt(np.sum(p0**2, axis=1))
	p0hat = (p0.T/p0abs).T
	dp = p-p0
	dp_para = np.sum(p0hat*dp, axis=1) #  dpz-1D array
	dp_perp2 = np.sum((dp - (p0hat.T*dp_para).T)**2, axis=1) # dpT2-1D array
	return dp_para.mean(), dp_perp2.mean()


def get_qhat(p, T):
	dt = 0.05
	qhatT3c = []
	qhatT3b = []
	e1 = event.event(medium={'type':'static', 'static_dt':dt},
				 LBT={'mu': 0.6}
				)
	e1.initialize_HQ(N_charm=10000, N_bottom=10000,
				init_flags={'type':'probe', 'E0':10.})
	for i in range(100):
		#t.append(e1.sys_time())
		print('{:1.2f}'.format(e1.sys_time()), '[fm/c]')
		e1.reset(pid=4, E0=(1.3**2+p**2)**0.5)
		e1.reset(pid=5, E0=(4.2**2+p**2)**0.5)
		pc0, xc0 = e1.HQ_hist(4)
		pb0, xb0 = e1.HQ_hist(5)
		status = e1.perform_hydro_step(StaticProperty={"Temp": T, "Vx":0.0, "Vy":0.0, "Vz":0.0})
		if e1.sys_time()>1.0:
			pc, xc = e1.HQ_hist(4)
			pb, xb = e1.HQ_hist(5)
			dpzc, dpT2c = decomp(pc0[:,1:], pc[:,1:])
			dpzb, dpT2b = decomp(pb0[:,1:], pb[:,1:])
			qhatT3c.append(dpT2c/dt/T**3/5.062)
			qhatT3b.append(dpT2b/dt/T**3/5.062)
		if not status:
			break
	return np.array([np.mean(qhatT3c), np.mean(qhatT3b)])

"""
with h5py.File("full-qhat.hdf5", 'a') as f:
	name = 'coll'
	if name in f:
		del f[name]
	gp = f.create_group(name)
	res = np.zeros([6,3,2])
	for iT,T in enumerate([0.16, 0.2, 0.3, 0.4, 0.5, 0.6]):
		for ip,p in enumerate([0.01, 5, 10]):
			res[iT, ip] = get_qhat(p, T)
	gp.attrs.create("T", [0.16, 0.2, 0.3, 0.4, 0.5, 0.6])
	gp.attrs.create("p", [0.01, 5, 10])
	gp.create_dataset('qhat', data=res)
"""
Ts = np.array([0.16, 0.2, 0.3, 0.4, 0.5, 0.6])
ps = np.array([0.01, 5, 10])
Ec = (1.3**2+ps**2)**0.5
Eb = (4.2**2+ps**2)**0.5
D = 0.4
x = 1.0
qhatT3_diff_c = D*(x+(1-x)/np.outer(Ts, Ec))
qhatT3_diff_b = D*(x+(1-x)/np.outer(Ts, Eb))
with h5py.File("full-qhat.hdf5", 'r') as f:
	for ip in range(3):
		plt.subplot(1,3,ip+1)
		plt.plot(Ts, 8*np.pi/(f['full/qhat'].value[:,:,0]+qhatT3_diff_c)[:, ip], 'c-')
		plt.plot(Ts, 8*np.pi/(f['coll/qhat'].value[:,:,0]+qhatT3_diff_c)[:, ip], 'c--')
		plt.plot(Ts, 8*np.pi/(f['full/qhat'].value[:,:,1]+qhatT3_diff_b)[:, ip], 'b-')
		plt.plot(Ts, 8*np.pi/(f['coll/qhat'].value[:,:,1]+qhatT3_diff_b)[:, ip], 'b--')
		plt.ylim(0,20)
	plt.show()
