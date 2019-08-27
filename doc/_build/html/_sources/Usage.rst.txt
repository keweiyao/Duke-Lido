Usage
====================================

There is a test.py file that demonstrates various setup of the code.
The code can work in either static medium or dynamical medium case.
The heavy quark evolution can be either a linear-boltzmann type evolution which uses LO-pQCD partonic cross-section or a Langevin type evolution.

First import the 'event' module

.. code-block:: python

  import event
  
  
To declare an event module, we need a bunch of options
They are split into two categories
Event options: 
1. Medium options: options are organized in a dictionary by { 'keyward' : value }
For example, a medium-option dictionary config a simulation with a static medium

.. code-block:: python
  
  static_config = { 'type'      : 'static',    
                    'static_dt' : 1.  }

The following medium information is to be used later

.. code-block:: python
  
  box_info = { 'Temp'  : 0.3, 
               'Vx'    : 0.0, 
               'Vy'    : 0.0, 
               'Vz'    : 0.0   }

Another example, a medium-option dictionary config a dynamical medium
with the path to the hydro-history file 

.. code-block:: python

  import sys
  hydro_history_filepath = sys.argv[1]
  dynamic_config = {  'type'      : 'dynamic', 
                      'hydrofile' : hydro_history_filepath    }

Physics option
1. An example for linear-Boltzmann evolution

.. code-block:: python
  
  LBT_config = {  'physics' : 'LBT',
                  '2->2'    : True,
                  '2->3'    : True,
                  '3->2'    : True,
                  'Nf'      : 3,
                  'mass'    : 1.3 }  
                                
2. An example for Langevin evolution

.. code-block:: python
  
  LGV_config = {  'physics'   : 'LGV',
                  'dt_lrf'    : 0.02,
                  'elastic'   : True,
                  'Einstein'  : True,
                  'Nf'        : 3,
                  'mass'      : 1.3 } 

Initialization option
Initizlize HQ in a static box :math:`[-L, L]^3`

.. code-block:: python

  import numpy
  box_init = {  'type'  : 'box',
                'L'     : 10.,
                'pmax'  : 10.   }
                
  TAA = numpy.loadtxt(sys.argv[2]).T
  
  realistic_init =  { 'type'          : 'A+B',
                      'sample power'  : 1.,
                      'pTmin'         : 0.1,
                      'pTmax'         : 70.,
                      'ymin'          : -1.,
                      'ymax'          : 1.,
                      'TAB'           : TAA,
                      'dxy'           : 0.1   }

.. code-block:: python
  
  e1 = event.event(   medium_flags=dynamic_config , 
                     physics_flags=LGV_config   )

  e1.initialize_HQ(   NQ=200000,
                      init_flags=realistic_init   )

Run Model

.. code-block:: python

  import h5py
  f = h5py.File("particle.hdf5", 'w')
  Init_pT = e1.Init_pT()
  f.create_dataset('init_pT', data=Init_pT)

  for i in range(500):
      print("t = %1.2f [fm/c]"%e1.sys_time() )
      status = e1.perform_hydro_step()#StaticPropertyDictionary=box_info)
      
      if i%5 == 0:
          dsp, dsx = e1.HQ_hist()
          f.create_dataset('p-%d'%i, data=dsp)
          
      if not status:
          break
          
  f.close()
