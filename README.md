The Duke Linearized-Boltzmann-and-Diffusion Partonic Transport Model 
====================================================================


#1 Install Pythia8
------------------

The collider mode study requires Pythia8. If you already have Pythia8 installed and have environment variables set in ``$HOME/.bashrc``, e.g,

```bash
   export PYTHIA8DATA=`pythia8-config --xmldoc`
   export PYTHIA8_DIR=`pythia8-config --prefix`
```

Otherwise please download [pythia8235](http://home.thep.lu.se/~torbjorn/pythiaaux/present.html) and install accordingly.


#2 Compile LIDO
-------------------

Requirements: ``c++11, gsl, hdf5, boost`` and make sure ``Pythia8`` is in your system path
Run cmake with the ``pythia8`` option turned on and then compile and install. Within the LIDO folder run the following:

```bash
   mkdir build && cd build
   cmake -Dpythia8=on ..
   make -j$(nproc)
   make install
```

This will install the package to the default localtion ``$HOME/.local/``. Or you can specifiy the install path by running ``cmake -DCMAKE_INSTALL_PREFIX:PATH=<your/install/path>`` instead.

#3 Download test hydro profiles.
----------------------------------

The LIDO transport model couples the jet partons to a hydrodynamic background. We have a few pregenerated hydrofiles for AuAu at 200 GeV. Download them by 

```bash
   ./get_test_hydro_profiles.sh
```

#4 Generate Tables
--------------------

LIDO will tabulate all those 2->2 and 2->3 and 1->2 cross-sections and rates. This can take sometime. But for every set of coupling/screening parameters, this only needs to be performed once. You can find the ``lido_settings.xml`` in the main folder which contains table information for all different scattering channles. To make table (assume you have install all LIDO excutables and add them to the system path)

```bash
    Lido-TabGen -s <path/to/lido_setting.xml> -t ./table.h5 --muT=1.5
```

will compute the tables with coupling scale parameters ``mu = 1.5*pi*T`` and output to the file ``table.h5``


