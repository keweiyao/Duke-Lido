The Duke Linearized-Boltzmann-and-Diffusion Partonic Transport Model  
====================================================================

Basic installation
------------------

Requirements: c++11, libraries: gsl, hdf5, and boost
Make and compile by running

```bash
   mkdir build && cd build
   cmake  ..
   make -j$(nproc) install
```

This will install the package to the default localtion ``$HOME/.local/``. Or you can specifiy the install path by running ``cmake -DCMAKE_INSTALL_PREFIX:PATH=<your/install/path>`` instead.

After this, you can also find an example under ``build/`` or ``<your/install/path/>``
For example, the ``analytic-medium`` executable run the model in a medium with a power-law temperature profile. ``./analytic-medium --help`` prints the options.

```bash
   ./analytic-medium -s <path-to-settings> -t <path-to-table> <more options> 
```

If the table does not exists, it will first take a few minutes to generate. The setting file is located in `examples/` or have already been installed to the install path under `share/`.

Install with Pythia8
--------------------

Some examples requires Pythia8 (prefer pythia8235). If you already have Pythia8 installed and have set the environment variables in your, e.g., ``$HOME/.bashrc``,

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

Coupled to hydrodynamic output
------------------------------

To be explained






