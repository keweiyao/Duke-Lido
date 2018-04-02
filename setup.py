#!/usr/bin/env python3
from setuptools import setup, Extension
from Cython.Build import cythonize
from glob import glob
import numpy as np
import os

# Remove the "-Wstrict-prototypes" compiler option, which isn't valid for C++.
import distutils.sysconfig
cfg_vars = distutils.sysconfig.get_config_vars()
for key, value in cfg_vars.items():
    if type(value) == str:
        cfg_vars[key] = value.replace("-Wstrict-prototypes", "")

libs=[path for path in os.environ['LD_LIBRARY_PATH'].split(':') if path]
print(os.environ.get('LD_LIBRARY_PATH'))
src = glob("./src/*.cpp")
#-------------HqEvo Module------------------
fileLBT = ['cython/event.pyx'] + src
modules = [
        Extension('event', 
        		 sources=fileLBT, 
        		 include_dirs=[np.get_include()],
        		 language="c++",
				 library_dirs=libs,
        		 extra_compile_args=["-std=c++11", '-fPIC'],
        		 libraries=["m", "gsl", "gslcblas", "boost_log", "boost_filesystem", "hdf5", "hdf5_cpp"])
]

setup(
        ext_modules=cythonize(modules),
	data_files=[('/share/event/', ['settings.xml'])]
)
