# cython: c_string_type=str, c_string_encoding=ascii, boundscheck=False
from libcpp cimport bool
from libcpp.map cimport map
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
from cython.parallel import parallel, prange

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
			print("works in static meidum mode!")
			print("Medium property can be specified step by step")
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

#--------alpha_s--------------------------------
cdef extern from "../src/matrix_elements.h":
	cdef double alpha_s(double Q2, double T)
falphas = open("alphas.dat", 'w')
#---------------------C++ lib for LBT and diffusion----------------
cdef extern from "../src/workflow.h":
	cdef struct particle:
		int pid;
		bool freezeout;
		double mass;
		fourvec x;
		fourvec p;
		bool has_k_rad, has_k_abs;
		double t_rad, t_abs;
		fourvec k_rad, k_abs;
		fourvec p0;
		vector[double] vcell;
		void freestream(double dt);
		double Tf;
		int resum_counts;

	cdef void initialize(string mode, string path, double mu, double alphafix)
	cdef int update_particle_momentum_Lido(double dt, double temp,
				vector[double] v3cell, particle)
	cdef int update_particle_momentum_HT(double dt, double temp,
				vector[double] v3cell, particle)
	cdef int gluon_elastic_scattering(double dt, double temp, vector[double] v3cell, fourvec incomping_p, fourvec & outgoing_p);
	cdef vector[double] probe_test(double M, double E0, double T, double dt, int Nsteps,
				int Nparticles, string mode, double mu, double alphafix);

	cdef vector[double] Bjorken_test(double E0, double T0, double t0, double dt, int Nsteps, int Nparticles, string mode, double mu, double const_alphas);

def probe_run(M, E0, T, dt=0.05, Nsteps=100, Nparticles=10000, mode="old", mu=1.0, alphafix=-1.0):
	dE = probe_test(M, E0, T, dt, Nsteps, Nparticles, mode, mu, alphafix)
	return dE;

def Bjorken_run(E0, T, t0=0.6, dt=0.05, Nsteps=100, Nparticles=10000, mode="old", mu=1.0, alphafix=-1.0):
	dE = Bjorken_test(E0, T, t0, dt, Nsteps, Nparticles, mode, mu, alphafix)
	return dE;

cdef extern from "../src/Langevin.h":
	cdef void initialize_transport_coeff(double A, double B)
	cdef void Ito_update(double dt, double M, double temp, vector[double] v3cell,
						const fourvec & pIn, fourvec & pOut)

cdef vector[double] regulate_v(vector[double] v):
	vabs2 = v[0]**2 + v[1]**2 + v[2]**2
	if vabs2 >= little_below_one**2:
		scale = 1.0/sqrt(vabs2)/little_above_one
		v[0] *= scale
		v[1] *= scale
		v[2] *= scale
	return v

cdef class event:
	cdef object hydro_reader, fs_reader
	cdef map[int, vector[particle]] HQ_list
	cdef str mode, transport
	cdef double Tc
	cdef double tau0, tau
	cdef bool lgv

	def __cinit__(self, preeq=None, medium=None,
			LBT=None, LGV=None, Tc=0.154):
		self.mode = medium['type']
		self.hydro_reader = Medium(medium_flags=medium)
		self.tau0 = self.hydro_reader.init_tau()
		self.HQ_list[4] = vector[particle]()
		self.HQ_list[5] = vector[particle]()
		if preeq is not None:
			self.fs_reader = Medium(medium_flags=preeq)
			self.tau0 = self.fs_reader.init_tau()
		self.tau = self.tau0
		self.lgv = False
		self.Tc = Tc

		# initialize LBT
		setting_path = os.environ['XDG_DATA_HOME']+"/event/settings.xml"
		print(setting_path)
		if not os.path.exists("table.h5"):
			initialize("new", setting_path, LBT['mu'], -1.0)
		else:
			initialize("old", setting_path, LBT['mu'], -1.0)

		# initialize LGV
		if LGV is not None:
			if LGV['A']>1e-9 and LGV['B'] > 1e-9:
				initialize_transport_coeff(LGV['A'], LGV['B'])
				self.lgv = True


	# The current time of the evolution.
	def sys_time(self) :
		return self.tau

	# Initilization
	cpdef initialize_HQ(self, N_charm, N_bottom, init_flags):
		cdef double x, y, z, s1, s2, mass
		cdef int pid
		# for A+B:
		#cdef double pT, phipt, rapidity, mT, t0
		cdef double Emax, ymax, pTmax, pTmin
		# for box:
		cdef double p, cospz, sinpz
		cdef double pmax, L
		cdef vector[particle].iterator it
		for pid, mass in zip([4,5],[1.3, 4.2]):
			self.HQ_list[pid].clear()
			NQ = N_charm if pid == 4 else N_bottom
			self.HQ_list[pid].resize(NQ) # NQ charm quark and NQ bottom quark, we don't need so many bottom quark
			if init_flags['type'] == 'A+B':
				print("Initialize for dynamic medium")
				HQ_xy_sampler = XY_sampler(init_flags['TAB'],
							   init_flags['dxy'],
							   init_flags['b'])
				logpTmax = log(init_flags['pTmax'])
				logpTmin = log(init_flags['pTmin'])
				Emax = init_flags['Emax']
				
				print("Heavy quarks are freestreamed to {} fm/c".format(self.tau0))
				it = self.HQ_list[pid].begin()
				X = []
				Y = []
				while it != self.HQ_list[pid].end():
					# Uniformly sample log(pT), phi, and ny = y/ymax
					# ymin, ymax are determined by the max-mT
					pT = np.exp(np.random.uniform(logpTmin, logpTmax))
					mT = sqrt(pT**2 + mass**2)
					phipt = np.random.uniform(0, 2.*np.pi)
					ymax = np.min([np.arccosh(Emax/mT), 3.0])
					rapidity = np.random.uniform(-ymax, ymax)
					pcharm = [mT*cosh(rapidity), pT*cos(phipt), \
							   pT*sin(phipt), mT*sinh(rapidity)]
					# sample initial x-y from TRENTo at time = 0+
					x, y, s1, s2 = HQ_xy_sampler.sample_xy()
					r0 = [0.0, x*fmc_to_GeV_m1, y*fmc_to_GeV_m1, 0.0]
					t0 = self.tau0/sqrt(1. - (pcharm[3]/pcharm[0])**2)*fmc_to_GeV_m1
					X.append(x)
					Y.append(y)
					# Initialize positional space at tau = 0+
					for i in range(4):
						deref(it).p.a[i] = pcharm[i]
						deref(it).p0.a[i] = pcharm[i]
						deref(it).x.a[i] = r0[i]
					deref(it).mass = mass
					# free streaming to hydro starting time tau = tau0
					deref(it).freestream(t0)
					# set last interaction vertex (assumed to be hydro start time)
					deref(it).t_rad = t0
					deref(it).t_abs = t0
					# initialize others
					deref(it).freezeout = False
					deref(it).vcell = [0., 0., 0.]
					deref(it).Tf = 0.
					deref(it).pid = pid
					inc(it)
				# check the variance of the sampling
				stdx, stdy = np.std(X), np.std(Y)
				print("std(x,y) = {:1.3f}, {:1.3f} [fm]".format(stdx, stdy) )
			elif init_flags['type'] == 'probe':
				print("Initialize for probe test")
				E0 = init_flags['E0']
				it = self.HQ_list[pid].begin()
				p0 = [E0, 0, 0, sqrt(E0*E0-mass*mass)]
				r0 = [0.0, 0.0, 0.0, 0.0]
				while it != self.HQ_list[pid].end():
					for i in range(4):
						deref(it).p.a[i] = p0[i]
						deref(it).p0.a[i] = p0[i]
						deref(it).x.a[i] = r0[i]
					deref(it).mass = mass
					deref(it).t_rad = r0[0]
					deref(it).t_abs = r0[0]
					deref(it).freezeout = False
					deref(it).vcell = [0., 0., 0.]
					deref(it).Tf = 0.
					deref(it).pid = pid
					inc(it)
			elif init_flags['type'] == 'Box':
				print("Initialize for probe test")
				pmax = init_flags['pmax']
				it = self.HQ_list[pid].begin()
				r0 = [0.0, 0.0, 0.0, 0.0]
				while it != self.HQ_list[pid].end():
					pT = np.random.rand()*pmax
					phi = np.random.rand()*2*np.pi
					cosz = np.random.rand()*2 - 1.
					sinz = np.sqrt(1. - cosz**2)
					E = np.sqrt(mass**2 + pT**2)
					p0 = [E, pT*sinz*np.cos(phi), pT*sinz*np.sin(phi), pT*cosz]
					r0 = [0,0,0,0]
					for i in range(4):
						deref(it).p.a[i] = p0[i]
						deref(it).p0.a[i] = p0[i]
						deref(it).x.a[i] = r0[i]
					deref(it).mass = mass
					deref(it).t_rad = r0[0]
					deref(it).t_abs = r0[0]
					deref(it).freezeout = False
					deref(it).vcell = [0., 0., 0.]
					deref(it).Tf = 0.
					deref(it).pid = pid
					inc(it)
			else:
				raise ValueError("Initilaiztion mode not defined")
				exit()

	cpdef bool perform_fs_step(self):
		PyErr_CheckSignals()
		if self.mode != 'dynamic':
			raise Warning("Only dyanmic mode has a freestream step")
			return False
		self.fs_reader.load_next()
		status = self.fs_reader.hydro_status()

		self.tau += self.fs_reader.dtau()
		cdef double T, tau_now, smaller_dtau, T1,T2
		cdef vector[double] vcell
		vcell.resize(3)
		cdef vector[particle].iterator it
		for pid in [4,5]:
			it = self.HQ_list[pid].begin()
			while it != self.HQ_list[pid].end():
				# use smaller time step than hydro
				for substeps in range(4):
					smaller_dtau = self.fs_reader.dtau()/4.
					# only update HQ that are not freezeout yet
					if not deref(it).freezeout:
						###############################################################
						###############################################################
						# Get the cell temperature and velocity for this heavy quark, #
						###############################################################
						tau_now = sqrt(deref(it).x.t()**2 - deref(it).x.z()**2)/fmc_to_GeV_m1
						T, vcell[0], vcell[1], vcell[2] = \
								self.hydro_reader.interpF(tau_now,
								[deref(it).x.t()/fmc_to_GeV_m1, 
								 deref(it).x.x()/fmc_to_GeV_m1,
								 deref(it).x.y()/fmc_to_GeV_m1,
								 deref(it).x.z()/fmc_to_GeV_m1],
								['Temp', 'Vx', 'Vy', 'Vz'])
						#	ensure |v| < 1.
						vcell = regulate_v(vcell)
						self.perform_HQ_step(it, tau_now, smaller_dtau, T, vcell)
				inc(it)
		return status

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
		#update system clock
		self.tau += self.hydro_reader.dtau()

		cdef double T, tau_now, smaller_dtau
		cdef vector[double] vcell
		vcell.resize(3)
		cdef vector[particle].iterator it
		for pid in [4,5]:
			it = self.HQ_list[pid].begin()
			while it != self.HQ_list[pid].end():
				# use smaller time step than hydro
				for substeps in range(10):
					smaller_dtau = self.hydro_reader.dtau()/10.
					# only update HQ that are not freezeout yet
					if not deref(it).freezeout:
						###############################################################
						###############################################################
						# Get the cell temperature and velocity for this heavy quark, #
						###############################################################
						if self.mode == "dynamic":
							tau_now = sqrt(deref(it).x.t()**2 - deref(it).x.z()**2)/fmc_to_GeV_m1
						else:
							tau_now = deref(it).x.t()
						T, vcell[0], vcell[1], vcell[2] = \
								self.hydro_reader.interpF(tau_now,
								[deref(it).x.t()/fmc_to_GeV_m1, 
								 deref(it).x.x()/fmc_to_GeV_m1,
								 deref(it).x.y()/fmc_to_GeV_m1,
								 deref(it).x.z()/fmc_to_GeV_m1],
								['Temp', 'Vx', 'Vy', 'Vz'])
						#	ensure |v| < 1.
						vcell = regulate_v(vcell)
						self.perform_HQ_step(it, tau_now, smaller_dtau, T, vcell)
				inc(it)
		return status

	cdef perform_HQ_step(self, vector[particle].iterator it,
						double tau_now, double dtau,
						double T, vector[double] vcell):
		PyErr_CheckSignals()
		cdef vector[fourvec] final_state
		cdef int channel

		###############################################################
		###############################################################
		# 	if T<Tc, label it as "freezout=True"			          #
		###############################################################
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

		deref(it).freestream(dt_lab*fmc_to_GeV_m1)
		channel = update_particle_momentum_Lido(
				dt_lab*fmc_to_GeV_m1, # evolve for this time
				T, vcell, deref(it))

		## if langevin is on, additional modification to momentum after LBT
		cdef fourvec pOut
		if self.lgv:
			Ito_update(dt_lab*fmc_to_GeV_m1, deref(it).mass, T, vcell,
						deref(it).p, pOut)
			deref(it).p = pOut
 

	cpdef HQ_hist(self, pid):
		cdef vector[particle].iterator it = self.HQ_list[pid].begin()
		cdef vector[ vector[double] ] p, x
		p.clear()
		x.clear()
		cdef fourvec ix, ip
		while it != self.HQ_list[pid].end():
			ip = deref(it).p
			ix = deref(it).x
			p.push_back([ip.t(),ip.x(),ip.y(),ip.z()])
			x.push_back([ix.t()/fmc_to_GeV_m1,ix.x()/fmc_to_GeV_m1,ix.y()/fmc_to_GeV_m1,ix.z()/fmc_to_GeV_m1])
			inc(it)
		return np.array(p), np.array(x)

	cpdef reset(self, int pid, double E0=10.):
		cdef vector[particle].iterator it = self.HQ_list[pid].begin()
		cdef double p0, rescale
		while it != self.HQ_list[pid].end():
			p0 = sqrt(E0**2 - deref(it).mass**2)
			rescale = p0/sqrt(deref(it).p.x()**2 + deref(it).p.y()**2 + deref(it).p.z()**2 )
			deref(it).p.a[1] = deref(it).p.x()*rescale
			deref(it).p.a[2] = deref(it).p.y()*rescale
			deref(it).p.a[3] = deref(it).p.z()*rescale
			deref(it).p.a[0] = E0
			inc(it)

	cpdef output_oscar(self, pid, filename):
		cdef vector[particle].iterator it = self.HQ_list[pid].begin()
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
				eventhead.write([1, self.HQ_list[pid].size(), 0.001, 0.001, 1, 1, 1])\
				+'\n')
			while it != self.HQ_list[pid].end():
				f.write(line.write([i, deref(it).pid,
					deref(it).p.x(),deref(it).p.y(),
					deref(it).p.z(),deref(it).p.t(),
					deref(it).mass,
					deref(it).x.x()/fmc_to_GeV_m1,deref(it).x.y()/fmc_to_GeV_m1,
					deref(it).x.z()/fmc_to_GeV_m1,deref(it).x.t()/fmc_to_GeV_m1,
					deref(it).Tf,
					deref(it).vcell[0], deref(it).vcell[1], deref(it).vcell[2],
					deref(it).p0.x(), deref(it).p0.y(),
					deref(it).p0.z(), deref(it).p0.t(),
					0., 0.])+'\n')
				i += 1
				inc(it)
