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

After this, you can also find an example under build/ or you/install/path/
```bash
   cp path/to/settings.xml ./
   ./example1 new
```

It should generate elastic scattering table and then calculate an energy loss for you.


