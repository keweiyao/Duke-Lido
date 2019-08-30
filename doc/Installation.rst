Installation
=================

To install, download the code from github and create a build directory

.. code::

  git clone -b https://github.com/keweiyao/Duke-Lido.git
  cd Duke-Lido
  mkdir build && cd build


The requirements are : a :code:`c++11` support compiler, :code:`cmake`>=3.4, c++ :code:`boost` library > 1.50, :code:`gsl` library, and c++ :code:`hdf5` library.
The default cmake also builds with :code:`Pythia8` support, but you can turn-off :code:`Pythia8` linking and disable the dependent exectuable by,

.. code::

  cmake -Dpythia8=off ..

otherwise,

.. code::

  cmake ..  





