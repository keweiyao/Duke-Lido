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
	cdef struct pregluon:
		fourvec p0, k1, kn;
		double t0, T0;
		double local_mfp;
		bool is_vac;
	cdef struct particle:
		int pid;
		bool freezeout;
		double mass;
		fourvec x;
		fourvec p;
		vector[pregluon] radlist, abslist;
		fourvec p0;
		vector[double] vcell;
		double Tf;
		double weight;
		void freestream(double dt);

	cdef void initialize(string mode, string path, double mu, double alphafix, double A, double B)
	cdef int update_particle_momentum_Lido(double dt, double temp,
				vector[double] v3cell, particle)

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
	cdef vector[particle] HQ_list
	cdef str mode
	cdef double Tc
	cdef double tau0, tau

	def __cinit__(self, preeq, medium, LIDO, Tc=0.154):
		self.mode = medium['type']
		self.HQ_list.clear()
		
		self.hydro_reader = Medium(medium_flags=medium)
		self.tau0 = self.hydro_reader.init_tau()
		
		self.fs_reader = Medium(medium_flags=preeq)
		self.tau0 = self.fs_reader.init_tau()
		
		self.tau = self.tau0
		self.Tc = Tc

		# initialize
		setting_path = os.environ['XDG_DATA_HOME']+"/event/settings.xml"
		print(setting_path)
		table_mode = "new" if not os.path.exists("table.h5") else "old"
		initialize(table_mode, setting_path, 
			LIDO['mu'], LIDO['afix'], LIDO['A'], LIDO['B'])


	# The current time of the evolution.
	def sys_time(self) :
		return self.tau

	# Initilization
	cpdef initialize_HQ(self, init_flags):
		cdef double x, y, z, s1, s2, tfs
		cdef particle H;
		cdef pregluon G;
		cdef vector[particle].iterator it
		self.HQ_list.clear()
		if init_flags['type'] == 'Pythia+Trento':
			print("Initialize four-mometa from a Pythia output")
			print("Initialize position from Trento T_AB")
			print("Heavy quarks are freestreamed to {} fm/c".format(self.tau0))

			HQ_xy_sampler = XY_sampler(init_flags['TAB'],
						   init_flags['dxy'],
						   init_flags['b'])
			
			Pid, M, E, Px, Py, Pz, Weight = np.loadtxt(init_flags['PythiaInput']).T
			for i, (pid, m, e, px, py, pz, w) in enumerate(zip(Pid, M, E, Px, Py, Pz, Weight)):
				if pid == 4:
					H.radlist.clear(); H.abslist.clear();	
					# sample initial x-y from TRENTo at time = 0+
					x, y, s1, s2 = HQ_xy_sampler.sample_xy()
					tfs = self.tau0/sqrt(1. - (pz/e)**2 + 1e-15)*fmc_to_GeV_m1
					# Initialize positional space at tau = 0+
					H.pid = pid
					H.mass = m			
					H.x.a[0] = 0.0; H.x.a[1] = x*fmc_to_GeV_m1;
					H.x.a[2] = y*fmc_to_GeV_m1; H.x.a[3] = 0.0;
					H.p.a[0] = e; H.p.a[1] = px;
					H.p.a[2] = py; H.p.a[3] = pz;
					# free streaming to hydro starting time tau = tau0
					H.freestream(tfs)
					# initialize others
					H.freezeout = False
					H.vcell = [0., 0., 0.]
					H.Tf = 0.
					H.weight = w
					self.HQ_list.push_back(H)
				else:
					G.p0 = self.HQ_list.back().p
					G.k1.a[0] = e; G.k1.a[1] = px;
					G.k1.a[2] = py; G.k1.a[3] = pz;
					G.kn = G.k1
					G.t0 = 0.0
					G.T0 = 0.0
					G.local_mfp = 0.0
					G.is_vac = True
					self.HQ_list.back().radlist.push_back(G)
		else:
			print("Initialization type {:s} for HQ does not exist".format(init_flags['type']))
			exit()

	cpdef bool perform_fs_step(self):
		PyErr_CheckSignals()
		if self.mode != 'dynamic':
			raise Warning("Only dyanmic mode has a freestream step")
			return False
		self.fs_reader.load_next()
		status = self.fs_reader.hydro_status()
		self.tau += self.fs_reader.dtau()

		cdef double T, tau_now
		cdef vector[double] vcell
		vcell.resize(3)
		cdef vector[particle].iterator it = self.HQ_list.begin()
		while it != self.HQ_list.end():
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
						self.fs_reader.interpF(tau_now,
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

	cpdef bool perform_hydro_step(self, StaticProperty={}):
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
		it = self.HQ_list.begin()
		while it != self.HQ_list.end():
			# use smaller time step than hydro
			for substeps in range(4):
				smaller_dtau = self.hydro_reader.dtau()/4.
				# only update HQ that are not freezeout yet
				if not deref(it).freezeout:
					###############################################################
					###############################################################
					# Get the cell temperature and velocity for this heavy quark, #
					###############################################################
					if self.mode == "dynamic":
						tau_now = sqrt(deref(it).x.t()**2 - deref(it).x.z()**2)/fmc_to_GeV_m1
					else:
						tau_now = deref(it).x.t()/fmc_to_GeV_m1
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
			t_m_zvz = deref(it).x.t()/fmc_to_GeV_m1 - deref(it).x.z()/fmc_to_GeV_m1*vz
			one_m_vz2 = 1. - vz*vz
			dtau2 = dtau*(dtau+2*tau_now)

			dt_lab = \
				(sqrt(t_m_zvz**2+one_m_vz2*dtau2) - t_m_zvz)/one_m_vz2
		else:
			dt_lab = dtau

		###############################################################
		###############################################################
		# update heavy quark status,                                  #
		# and return the full final states                            #
		###############################################################
		# time should be in GeV^-1 in the update function !!!

		deref(it).freestream(dt_lab*fmc_to_GeV_m1)
		update_particle_momentum_Lido(
				dt_lab*fmc_to_GeV_m1, # evolve for this time, Langevin included
				T, vcell, deref(it))
 

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
			x.push_back([ix.t()/fmc_to_GeV_m1,ix.x()/fmc_to_GeV_m1,ix.y()/fmc_to_GeV_m1,ix.z()/fmc_to_GeV_m1])
			inc(it)
		return np.array(p), np.array(x)

	cpdef reset(self, double E0=10.):
		cdef vector[particle].iterator it = self.HQ_list.begin()
		cdef double p0, rescale
		while it != self.HQ_list.end():
			p0 = sqrt(E0**2 - deref(it).mass**2)
			rescale = p0/sqrt(deref(it).p.x()**2 + deref(it).p.y()**2 + deref(it).p.z()**2 )
			deref(it).p.a[1] = deref(it).p.x()*rescale
			deref(it).p.a[2] = deref(it).p.y()*rescale
			deref(it).p.a[3] = deref(it).p.z()*rescale
			deref(it).p.a[0] = E0
			inc(it)

	cpdef output_oscar(self, pid, filename):
		cdef vector[particle].iterator it
		cdef size_t i=0, Ntot=0
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

			it = self.HQ_list.begin()
			while it != self.HQ_list.end():
				if deref(it).pid == pid:
					Ntot += 1
				inc(it)
			f.write(
				eventhead.write([1, Ntot, 0.001, 0.001, 1, 1, 1])\
				+'\n')
			it = self.HQ_list.begin()
			while it != self.HQ_list.end():
				if deref(it).pid == pid:
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
						deref(it).weight, 0.])+'\n')
					i += 1
				inc(it)
