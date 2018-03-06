# cython: c_string_type=str, c_string_encoding=ascii
from libcpp cimport bool
from libcpp.vector cimport vector
from libcpp.string cimport string
from libc.stdlib cimport rand, RAND_MAX
from libc.math cimport *
from libc.stdlib cimport malloc, free
from cython.operator cimport dereference as deref, preincrement as inc
from cpython.exc cimport PyErr_CheckSignals
import numpy as np
cimport numpy as np
import h5py
import os

import fortranformat as ff

cdef double GeV_m1_to_fmc = 0.197
cdef double fmc_to_GeV_m1 = 5.026
cdef double little_below_one = 1. - 1e-6
cdef double little_above_one = 1. + 1e-6		

#-----------Hydro reader class--------------------------------------
cdef inline double finterp(double *** c, double rx, double ry, double rz):
	cdef double result = 0.0
	cdef size_t i, j, k
	cdef double vx[2]
	vx[0] = 1. - rx; vx[1] = rx
	cdef double vy[2]
	vy[0] = 1. - ry; vy[1] = ry
	cdef double vz[2]
	vz[0] = 1. - rz; vz[1] = rz
	for i in range(2):
		for j in range(2):
			for k in range(2):
				result += c[i][j][k]*vx[i]*vy[j]*vz[k]
	return result

cdef class Medium:
	cdef object _f, _step_keys, info_keys, static_property, _tabs
	cdef str _mode
	cdef public size_t _Nx, _Ny, _step_key_index
	cdef public double _dx, _dy, _xmin, _ymin, _xmax, _ymax, _tstart, _dt, _tnow
	cdef double T_static
	cdef double *** interpcube
	cdef bool status

	def __cinit__(self, medium_flags):
		self._mode = medium_flags['type']

		if self._mode == 'static':
			print "works in static meidum mode!"
			print "Medium property can be specified step by step"
			self._Nx = 0
			self._Ny = 0
			self._dx = 0.
			self._dy = 0.
			self._tstart = 0.0
			self._dt = medium_flags['static_dt']
			self._xmin = 0.
			self._xmax = 0.
			self._ymin = 0.
			self._ymax = 0.
			self._tnow = self._tstart - self._dt
			self.status = True
		elif self._mode == 'dynamic':
			hydrofilename = medium_flags['hydrofile']
			if hydrofilename == None:
				raise ValueError("Need hydro history file")
			else:
				self._f = h5py.File(hydrofilename, 'r')

			self._step_keys = list(self._f['Event'].keys())
			self._step_key_index = 0
			self.status = True
			self.info_keys = self._f['Event'][self._step_keys[0]].keys()
			self._Nx = self._f['Event'].attrs['XH'] - self._f['Event'].attrs['XL'] + 1
			self._Ny = self._f['Event'].attrs['YH'] - self._f['Event'].attrs['YL'] + 1
			self._dx = self._f['Event'].attrs['DX']
			self._dy = self._f['Event'].attrs['DY']
			self._tstart = self._f['Event'].attrs['Tau0']
			self._dt = self._f['Event'].attrs['dTau']
			self._xmin = self._f['Event'].attrs['XL']*self._dx
			self._xmax = self._f['Event'].attrs['XH']*self._dx
			self._ymin = self._f['Event'].attrs['YL']*self._dy
			self._ymax = self._f['Event'].attrs['YH']*self._dy
			self._tnow = self._tstart - self._dt

			self.interpcube = <double ***> malloc(2*sizeof(double**))
			for i in range(2):
				self.interpcube[i] = <double **> malloc(2*sizeof(double*))
				for j in range(2):
					self.interpcube[i][j] = <double *> malloc(2*sizeof(double))
		else:
			raise ValueError("Medium mode not implemented.")

	cpdef init_tau(self):
		return self._tstart
	cpdef hydro_status(self):
		return self.status
	cpdef dtau(self):
		return self._dt
	cpdef boundary(self):
		return self._xmin, self._xmax, self._ymin, self._ymax

	cdef frame_inc_unpack(self):
		self._tabs = {}
		cdef vector[vector[vector[double]]] buff
		buff.resize(2)
		for key2 in self.info_keys:
			for i in range(2):
				key1 = self._step_keys[self._step_key_index+i]
				buff[i] = self._f['Event'][key1][key2].value
			self._tabs.update({key2:buff})
		self._step_key_index += 1
		if self._step_key_index == len(self._step_keys) - 1:
			self.status = False

	cpdef load_next(self, StaticProperty=None):
		if self._mode == "dynamic":
			self.frame_inc_unpack()
		elif self._mode == 'static':
			if StaticProperty == None:
				raise ValueError("Requires static meidum properties")
			else:
				self.static_property = StaticProperty
		else:
			raise ValueError("Medium mode not implemented.")
		self._tnow += self._dt

	cpdef get_current_frame(self, key):
		return np.array(self._tabs[key][0])

	cpdef interpF(self, double tau, xvec, keys):
		cdef double rt, nx, ny, rx, ry, gamma, buff, vz
		cdef int ix, iy, i, j, k
		if self._mode == "static":
			return [self.static_property[key] for key in keys]
		if self._mode == "dynamic":
			if xvec[1] < self._xmin or xvec[1] > self._xmax \
				or xvec[2] < self._ymin or xvec[2] > self._ymax:
				return [0. for key in keys]
			else:
				result = []
				rt = (tau - self._tnow)/self._dt
				nx = (xvec[1] - self._xmin)/self._dx
				ny = (xvec[2] - self._ymin)/self._dy
				ix = <int>floor(nx)
				rx = nx - ix
				iy = <int>floor(ny)
				ry = ny - iy
				vz = xvec[3]/xvec[0]
				for key in keys:
					if key == 'Vz':
						result.append(vz)
					else:
						for k in range(2):
							for i in range(2):
								for j in range(2):
									self.interpcube[k][i][j] = \
										self._tabs[key][k][ix+i][iy+j]
						buff = finterp(self.interpcube, rt, rx, ry)
						if key == 'Vx' or key == 'Vy':
							gamma = 1.0/sqrt(1.0-vz*vz)
							buff /= gamma
						result.append(buff)
				return result

#-----------Production vertex sampler class-------------------------
cdef class XY_sampler:
	cdef np.ndarray Taa, IntTaa
	cdef double dxy, b
	cdef size_t Nx, Ny
	def __cinit__(self, Taa, dxy, b):
		self.b = b
		self.dxy = dxy
		self.Nx, self.Ny = Taa.shape
		self.Taa = Taa.reshape(-1)
		self.IntTaa = np.zeros_like(self.Taa, dtype=np.double)
		cdef double tot = self.Taa.sum()
		cdef int i
		cdef double dT
		for i, dT in enumerate(self.Taa):
			self.IntTaa[i] = self.IntTaa[i-1] + dT/tot
	cpdef sample_xy(self):
		cdef double r = np.random.rand()
		cdef int index = np.searchsorted(self.IntTaa, r)
		cdef double nx = np.floor((index-1.)/self.Ny)
		cdef double ny = index - 1 - nx*self.Ny
		cdef double x, y, s1, s2
		nx += np.random.rand()
		ny += np.random.rand()
		# to be examined
		x, y = (nx - self.Nx/2.)*self.dxy, (ny - self.Ny/2.)*self.dxy
		s1 = sqrt(y**2 + (x+self.b/2.)**2)
		s2 = sqrt(y**2 + (x-self.b/2.)**2)
		return x, y, s1, s2

#---------------Particle data------------------------------------
cdef extern from "../src/lorentz.h":
	cdef struct scalar:
		double s
	cdef struct fourvec:
		double a[4]
		double t()
		double x()
		double y()
		double z()
		fourvec boost_back(double vx, double vy, double vz)
		fourvec boost_to(double vx, double vy, double vz)
		fourvec rotate_back(const fourvec p)

#---------------------C++ lib for LBT and diffusion----------------
cdef extern from "../src/workflow.h":
	cdef struct particle:
		int pid
		bool freezeout
		double mass
		fourvec x
		fourvec p
		double t_rad, t_absorb
		fourvec p0
		vector[double] vcell
		double Tf
		void freestream(double dt)
	cdef void initialize(string s)
	cdef int update_particle_momentum(double dt, double temp, 
				vector[double] v3cell, 
				double D_formation_t, fourvec incoming_p, 
				vector[fourvec] & FS);

cdef double proper_time(vector[double] x):
	return sqrt(x[0]**2 - x[3]**2)

cdef class event:
	cdef object hydro_reader, lbt, lgv
	cdef str mode, transport
	cdef double M, Tc
	cdef vector[particle] HQ_list
	cdef double tau0, dtau, tau

	def __cinit__(self, medium={"type":"static", "static_dt":0.05}, 
					LBT=None, LGV=None, Tc=0.154, M=1.3):
		# read medium config
		self.mode = medium['type']
		self.hydro_reader = Medium(medium_flags=medium)
		self.tau0 = self.hydro_reader.init_tau()
		self.dtau = self.hydro_reader.dtau()
		self.tau = self.tau0

		# load transport
		self.Tc = Tc
		self.M = M

		# initialize LBT
		if not os.path.exists("table.h5"):
			initialize("new")
		else:
			initialize("old")	
	# The current time of the evolution.
	def sys_time(self) :
		return self.tau

	# Initilization
	cpdef initialize_HQ(self, NQ, 
			init_flags={"type":'box', "pmax":10., "L":10.}):
		cdef double x, y, z, s1, s2
		# for A+B:
		#cdef double pT, phipt, rapidity, mT, t0
		cdef double ymax, pTmax, pTmin
		# for box:
		cdef double p, cospz, sinpz
		cdef double pmax, L
		cdef vector[particle].iterator it

		if init_flags['type'] == 'A+B':
			self.HQ_list.clear()
			self.HQ_list.resize(NQ)
			print("Initialize for dynamic medium")
			HQ_xy_sampler = XY_sampler(init_flags['TAB'],
									   init_flags['dxy'],
									   init_flags['b'])
			pTmax = init_flags['pTmax']
			pTmin = init_flags['pTmin']
			ymax = init_flags['ymax']

			print("Heavy quarks are freestreamed to {} fm/c".format(self.tau0))
			it = self.HQ_list.begin()
			X = []
			Y = []
			while it != self.HQ_list.end():
				# Uniformly sample pT, phi, and y
				pT = np.random.uniform(pTmin, pTmax)
				mT = sqrt(pT**2 + self.M**2)
				phipt = np.random.uniform(0, 2.*np.pi)
				rapidity = np.random.uniform(-ymax, ymax)
				pcharm = [mT*cosh(rapidity), pT*cos(phipt), \
							   pT*sin(phipt), mT*sinh(rapidity)]
				# sample initial x-y from TRENTo at time = 0+
				x, y, s1, s2 = HQ_xy_sampler.sample_xy()
				r0 = [0.0, x, y, 0.0]
				t0 = self.tau0/sqrt(1. - (pcharm[3]/pcharm[0])**2)
				X.append(x)
				Y.append(y)
				# charm:
				# Initialize positional space at tau = 0+
				for i in range(4):
					deref(it).p.a[i] = pcharm[i]
					deref(it).x.a[i] = r0[i]
				deref(it).mass = self.M
				# free streaming to hydro starting time tau = tau0
				deref(it).freestream(t0)
				# set last interaction vertex (assumed to be hydro start time)
				deref(it).t_rad = t0
				deref(it).t_absorb = t0
				# initialize others
				deref(it).freezeout = False
				deref(it).p0 = deref(it).p
				deref(it).vcell = [0., 0., 0.]
				deref(it).Tf = 0.
				deref(it).pid = 4
				inc(it)
			# check the variance of the sampling
			stdx, stdy = np.std(X), np.std(Y)
			print("std(x,y) = {:1.3f}, {:1.3f} [fm]".format(stdx, stdy) )

		elif init_flags['type'] == 'box':
			self.HQ_list.clear()
			self.HQ_list.resize(NQ)
			print "Initialize for box simulation"
			pmax = init_flags['pmax']
			L = init_flags['L']
			it = self.HQ_list.begin()

			while it != self.HQ_list.end():
				p = np.random.uniform(0, pmax)
				phipt = np.random.uniform(0, 2.*np.pi)
				cospz = little_below_one*np.random.uniform(-1., 1.)
				sinpz = sqrt(1.-cospz**2)
				p0 = [sqrt(p**2+self.M**2), p*sinpz*cos(phipt), \
								p*sinpz*sin(phipt), p*cospz]
				r0 = [0.0, x, y, z]
				for i in range(4):
					deref(it).p.a[i] = p0[i]
					deref(it).x.a[i] = r0[i]
				deref(it).t_rad = 0.0
				deref(it).t_absorb = 0.0
				deref(it).freezeout = False
				deref(it).p0 = deref(it).p
				deref(it).vcell = [0., 0., 0.]
				deref(it).Tf = 0.
				deref(it).pid = 4
				deref(it).mass = self.M
				x,y,z = np.random.uniform(-L,L,3)
				
				inc(it)
		else:
			raise ValueError("Initilaiztion mode not defined")
			exit()

	cpdef bool perform_hydro_step(self, 
			StaticProperty={"Temp": 0.4, "Vx":0.0, "Vy":0.0, "Vz":0.0}	):
		PyErr_CheckSignals()
		if self.mode == 'dynamic':
			self.hydro_reader.load_next()
		elif self.mode == 'static':
			self.hydro_reader.load_next(StaticProperty=StaticProperty)
		else:
			raise ValueError("medium mode not defined")
		status = self.hydro_reader.hydro_status()

		self.tau += self.dtau
		cdef double t, x, y, z, tauQ
		cdef vector[particle].iterator it
		it = self.HQ_list.begin()
		while it != self.HQ_list.end():
			if not deref(it).freezeout:
				for substeps in range(2):
					self.perform_HQ_step(it, self.dtau/2.)
			inc(it)
		return status

	cdef perform_HQ_step(self, vector[particle].iterator it, double dtau):
		PyErr_CheckSignals()
		cdef double t, tau_now, T, vabs2, scale
		cdef vector[double] vcell
		cdef vector[fourvec] final_state
		cdef int channel

		###############################################################
		###############################################################
		# Get the cell temperature and velocity for this heavy quark, #
		# ensure |v| < 1. if T<Tc, label it as "freezout=True"        #
		###############################################################
		if self.mode == "dynamic":
			tau_now = sqrt(deref(it).x.t()**2 - deref(it).x.z()**2)
		else:
			tau_now = deref(it).x.t()
		vcell.resize(3)
		T, vcell[0], vcell[1], vcell[2] = \
			self.hydro_reader.interpF(tau_now, 
				[deref(it).x.t(),deref(it).x.x(),deref(it).x.y(),deref(it).x.z()],
				['Temp', 'Vx', 'Vy', 'Vz'])
		vabs2 = vcell[0]**2 + vcell[1]**2 + vcell[2]**2
		if vabs2 >= little_below_one**2:
			scale = 1.0/sqrt(vabs2)/little_above_one
			vcell[0] *= scale
			vcell[1] *= scale
			vcell[2] *= scale
		if T <= self.Tc:
			deref(it).freezeout = True
			deref(it).Tf = T
			deref(it).vcell = vcell
			return

		###############################################################
		###############################################################
		# how long does it need to evolve to reach the next proper    #
		# time step                                                   #
		###############################################################
		cdef double vz, t_m_zvz, one_m_vz2, dtau2, dt_lab
		if self.mode == "dynamic":
			vz = deref(it).p.z()/deref(it).p.t()
			t_m_zvz = deref(it).x.t() - deref(it).x.z()*vz
			one_m_vz2 = 1. - vz*vz
			dtau2 = dtau*(dtau+2*tau_now)
		 
			dt_lab = \
				(sqrt(t_m_zvz**2+one_m_vz2*dtau2) - t_m_zvz)/one_m_vz2
		else:
			dt_lab = dtau
		###############################################################
		###############################################################
		# update heavy quark status, returns which scatterig channel, #
		# and return the full final states                            #
		###############################################################
		# time should be in GeV^-1 in the update function !!!
		channel = update_particle_momentum(
				dt_lab*fmc_to_GeV_m1, # evolve for this time 
				T, vcell, 	# fluid info
				(deref(it).x.t() - deref(it).t_rad)*fmc_to_GeV_m1, # LPM effect
				deref(it).p, # Initial probe momentum
				final_state # Final states
			)

		if channel < 0: # Nothing happens
			deref(it).freestream(dt_lab)
		elif channel == 0 or channel == 1: #elastic
			deref(it).freestream(dt_lab)
			deref(it).p = final_state[0]
		elif channel == 2 or channel == 3:
			deref(it).freestream(dt_lab)
			deref(it).p = final_state[0]
			deref(it).t_rad = deref(it).x.t()
		elif channel == 4 or channel == 5:
			deref(it).freestream(dt_lab)
			deref(it).p = final_state[0]
			deref(it).t_absorb = deref(it).x.t()
		else:
			raise ValueError("Unknown channel")

	cpdef HQ_hist(self):
		cdef vector[particle].iterator it = self.HQ_list.begin()
		cdef vector[ vector[double] ] p, x
		p.clear()
		x.clear()
		cdef fourvec ix, ip
		while it != self.HQ_list.end():
			ip = deref(it).p
			ix = deref(it).x
			p.push_back([ip.t(),ip.x(),ip.y(),ip.z()])
			x.push_back([ix.t(),ix.x(),ix.y(),ix.z()])
			inc(it)
		return np.array(p), np.array(x)

	cpdef output_oscar(self, filename):
		cdef vector[particle].iterator it = self.HQ_list.begin()
		cdef size_t i=0
		with open(filename, 'w') as f:
			head3 = ff.FortranRecordWriter(
                 '2(a8,2x),1x,i3,1x,i6,3x,i3,1x,i6,3x,a4,2x,e10.4,2x,i8')
			f.write('OSC1997A\n')
			f.write('final_id_p_x\n')
			f.write(head3.write(['lbt', '1.0alpha', 208, 82, \
								 208, 82, 'aacm', 1380, 1])+'\n')
			line = ff.FortranRecordWriter('i10, 2x, i10,19(2x,d12.6)')
			eventhead =ff.FortranRecordWriter(
					'i10,2x,i10,2x,f8.3,2x,f8.3,2x,i4,2x,i4,2X,i7')
			f.write(
				eventhead.write([1, self.HQ_list.size(), 0.001, 0.001, 1, 1, 1])\
				+'\n')
			while it != self.HQ_list.end():
				f.write(line.write([i, deref(it).pid,
					deref(it).p.x(),deref(it).p.y(),
					deref(it).p.z(),deref(it).p.t(),
					deref(it).mass,
					deref(it).x.x(),deref(it).x.y(),
					deref(it).x.z(),deref(it).x.t(),
					deref(it).Tf,
					deref(it).vcell[0], deref(it).vcell[1], deref(it).vcell[2],
					deref(it).p0.x(), deref(it).p0.y(),
					deref(it).p0.z(), deref(it).p0.t(),
					0., 0.])+'\n')
				i += 1
				inc(it)

