The Duke Linearized-Boltzmann-and-Diffusion Partonic Transport Model 
====================================================================


1. Install Pythia8
--------------------

Some examples require Pythia8. If you already have Pythia8 installed and have environment variables set in ``$HOME/.bashrc``, e.g,

```bash
   export PYTHIA8DATA=`pythia8-config --xmldoc`
   export PYTHIA8_DIR=`pythia8-config --prefix`
```

Then run cmake with the ``pythia8`` option turned on and then compile and install

```bash
   mkdir build && cd build
   cmake -Dpythia8=on ..
   make -j$(nproc) install
```

Otherwise please download [pythia8235](http://home.thep.lu.se/~torbjorn/pythiaaux/present.html) and install accordingly.

2. Install LIDO
------------------

Requirements: c++11, libraries: gsl, hdf5, and boost and make sure Pythia8 is in your system path
Make and compile by running

```bash
   mkdir build && cd build
   cmake  ..
   make -j$(nproc) install
```

This will install the package to the default localtion ``$HOME/.local/``. Or you can specifiy the install path by running ``cmake -DCMAKE_INSTALL_PREFIX:PATH=<your/install/path>`` instead.

