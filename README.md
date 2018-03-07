Scattering 
==========
-------------------------------------------------
New design for the Duke linear Boltzmann project.
-------------------------------------------------

Requirements: c++11, libraries: gsl, hdf5, and boost

To compile and install exe:

```bash
   mkdir build && cd build
   cmake ..
   make
   make install
```


To setup python interface, have Cython installed and then

```bash
   python3 setup --build_ext -b <local folder>
```
or 

```bash
   python3 setup install --prefix=<path to install>
```
