from setuptools import setup, Extension
from Cython.Build import cythonize
from glob import glob
import os

# Remove the "-Wstrict-prototypes" compiler option, which isn't valid for C++.
import distutils.sysconfig
cfg_vars = distutils.sysconfig.get_config_vars()
for key, value in cfg_vars.items():
    if type(value) == str:
        cfg_vars[key] = value.replace("-Wstrict-prototypes", "")

#-------------HqEvo Module------------------
fileLBT = [	'cython/pyScattering.pyx', 
			'src/matrix_elements.cpp']
modules = [
        Extension('pyScattering', 
        		 sources=fileLBT, 
        		 language="c++", 
        		 extra_compile_args=["-std=c++11", '-march=native', '-fPIC'],
        		 libraries=["m", "gsl", "gslcblas", "boost_filesystem", "hdf5", "hdf5_cpp"])
]


setup(
        ext_modules=cythonize(modules),
)
