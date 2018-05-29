#!/usr/bin/env python3
from scipy.spatial import cKDTree
import numpy as np
import matplotlib.pyplot as plt

ptype = [('x', np.float, 3), ('v', np.float, 3), ('M', np.float)]
def vrel(p1,p2):
	return p1['v']-p2['v']
def vcom_decompose(p1,p2):
	vcom = (p1['v']*p1['M']+p2['v']*p2['M'])/(p1['M']+p2['M'])
	v1com = p1['v']-vcom
	v2com = p2['v']-vcom
	return vcom, v1com, v2com
class solver:
	def __init__(self, N, x0, v0, sigma0, sigma):
		# initialize a almond shaped spatial distribution and isotropic velocity
		self.plist = np.zeros(N, dtype=ptype)
		for i, p in enumerate(self.plist):
			self.plist[i]['x'] = np.random.normal(0,1,3)*x0*np.array([1,0.2,1])
			self.plist[i]['v'] = v0*np.random.normal(0,1,3)
			self.plist[i]['M'] = 1.
		self.sigma0 = sigma0 # cross-section
		self.N = N
		self.update_tree()
	def hist(self,key):
		return np.array([p[key] for p in self.plist])
	def update_tree(self): 
		# Construct KD-tree to find nearist neighbor
		self.tree = cKDTree([p['x'] for p in self.plist])
		# update dimensionaless parameter, and dt
		stdx, stdy, stdz = np.std([p['x'] for p in self.plist], axis=0)
		stdvx, stdvy, stdvz = np.std([p['v'] for p in self.plist], axis=0)
		self.Leff = np.sqrt(stdx**2+stdy**2+stdz**2)
		self.Veff = np.sqrt(stdvx**2+stdvy**2+stdvz**2)
		self.mfp = (self.Leff)**3/self.N/self.sigma0
		self.r0 = self.Leff/20.
		self.dV = 4.*np.pi/3.*self.r0**3
		self.dt = np.min([self.mfp/self.Veff, self.dV/self.sigma0/self.Veff])/10.
	def print_characterist_number(self):
		print("mfp ~ ", self.mfp)
		print("Knudson number ~ ", self.mfp/self.Leff)
		print("dt = ", self.dt)		
	def evolve(self):
		self.update_tree()
		pairs = self.tree.query_pairs(self.r0)
		for (i1, i2) in pairs:
			p1, p2 = self.plist[i1], self.plist[i2]
			prob = np.sqrt((vrel(p1, p2)**2).sum())*self.dt*self.sigma0/self.dV
			if prob > 1.:
				print('prob = ', prob, 
						' is too large to suppress multiple collision')
			if np.random.rand() < prob:
				vcom, v1com, v2com = vcom_decompose(p1,p2)
				# isotropic scattering
				phi = np.random.rand()*np.pi*2
				cz = np.random.rand()*2.-1.
				sz = np.sqrt(1.-cz**2)
				cphi, sphi = np.cos(phi), np.sin(phi)
				direction = np.array([sz*cphi, sz*sphi, cz])
				self.plist[i1]['v'] = vcom + direction*np.sqrt((v1com**2).sum())
				self.plist[i2]['v'] = vcom - direction*np.sqrt((v2com**2).sum())
		# spatial transport
		for i, p in enumerate(self.plist):
			self.plist[i]['x'] += p['v']*self.dt

N = 10000
t = 0
s = solver(N, x0=.5, v0=1., sigma0=.001)
fig = plt.figure(figsize=(6,3))
for i in range(220):
	t += s.dt	
	if i%10 == 0:
		print(i)
		s.print_characterist_number()
		plt.clf()
		vx, vy, vz = s.hist('v').T
		rx, ry, rz = s.hist('x').T
		plt.subplot(1,2,1)
		H, a, b = np.histogram2d(rx, ry, bins=51, range=[[-5,5],[-5,5]])
		plt.contourf(H, extent=[-5,5,-5,5])
		plt.axis([-5,5,-5,5])
		plt.xlabel("x")
		plt.ylabel("y")
		plt.subplot(1,2,2)
		H, a, b = np.histogram2d(vx, vy, bins=51, range=[[-5,5],[-5,5]])
		plt.contourf(H, extent=[-5,5,-5,5])
		plt.axis([-5,5,-5,5])
		plt.xlabel("vx")
		plt.ylabel("vy")
		plt.subplots_adjust(bottom=0.15,left=0.15, right=0.95, wspace=0.25)
		plt.pause(.1)
	s.evolve()
plt.show()






