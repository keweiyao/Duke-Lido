# cython: c_string_type=str, c_string_encoding=ascii
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp cimport bool
from libc.stdlib cimport malloc, free
from libc.stdlib cimport rand, RAND_MAX
from libc.math cimport fmin
import cython
import os
import numpy as np
cimport numpy as np
import subprocess

cdef double GeV_to_Invfm = 5.068

cdef extern from "../src/matrix_elements.h":
	cdef void initialize_mD_and_scale(const unsigned int mDtype, const double scale)
	cdef double dX_Qq2Qq_dt(double t, void * params)
	cdef double dX_Qg2Qg_dt(double t, void * params)


cdef extern from "../src/lorentz.h":
	cdef struct scalar:
		double s
	cdef struct fourvec:
		double a[4]
		fourvec boost_back(double vx, double vy, double vz)
		fourvec boost_to(double vx, double vy, double vz)
		fourvec rotate_back(const fourvec p)
	
cdef extern from "../src/workflow.h":
	cdef struct particle:
		int pid
		fourvec x, p
		double trad
		void freestream(double dt)
	cdef void initialize(string s)
	cdef int update_particle_momentum(double dt, double temp, 
				vector[double] v3cell, 
				double D_formation_t, fourvec incoming_p, vector[fourvec] & FS);

	cdef void probe_test(double E0, double T, double dt, int Nsteps, 
				int Nparticles, string mode);	


initialize_mD_and_scale(0, 2.0)

def py_dX_Qq2Qq_dt(double t, np.ndarray[np.double_t, ndim=1] X):
	X = np.ascontiguousarray(X)
	return dX_Qq2Qq_dt(t, &X[0])

def py_dX_Qg2Qg_dt(double t, np.ndarray[np.double_t, ndim=1] X):
	X = np.ascontiguousarray(X)
	return dX_Qg2Qg_dt(t, &X[0])
	
def probe_run(E0, T, dt=0.05, Nsteps=100, Nparticles=10000, mode="old"):
	probe_test(E0, T, dt, Nsteps, Nparticles, mode)	
	

		



