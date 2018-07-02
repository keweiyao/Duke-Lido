Installation
=================

Requirements: c++11, libraries: gsl, hdf5, and boost

To compile and install exe:

.. code-block:: bash

  mkdir build && cd build
  cmake ..
  make
  make install


After this, you can also find an example under build/ or you/install/path/

.. code-block:: bash

  cp path/to/settings.xml ./
  ./example1 new


It should generate elastic scattering table and then calculate an energy loss for you.


To setup python interface, have Cython installed and then

.. code-block:: bash

  python3 setup --build_ext -b <local folder>

or 

.. code-block:: bash

  python3 setup install --prefix=<path to install>

