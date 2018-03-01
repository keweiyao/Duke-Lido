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

cdef double GeV_to_Invfm = 5.068

cdef extern from "../src/matrix_elements.h":
	cdef void initialize_mD_and_scale(const unsigned int mDtype, const double scale)

	cdef double dX_Qq2Qq_dt(double t, void * params)
	cdef double dX_Qg2Qg_dt(double t, void * params)

cdef extren from "../src/Xsection.h":
	cdef Xsection[size_t N, typename F]
		Xsection()
initialize_mD_and_scale(0, 2.0)

def py_dX_Qq2Qq_dt(double t, np.ndarray[np.double_t, ndim=1] X):
	X = np.ascontiguousarray(X)
	return dX_Qq2Qq_dt(t, &X[0])

def py_dX_Qg2Qg_dt(double t, np.ndarray[np.double_t, ndim=1] X):
	X = np.ascontiguousarray(X)
	return dX_Qg2Qg_dt(t, &X[0])

