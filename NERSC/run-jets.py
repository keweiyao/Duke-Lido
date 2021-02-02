#!/usr/bin/env python3

import argparse
from contextlib import contextmanager
import datetime
from itertools import chain, groupby, repeat
import logging
import math
import os
import pickle
import signal
import subprocess
import sys
import tempfile
import numpy as np
import h5py
from scipy.interpolate import interp1d
from scipy.ndimage.measurements import center_of_mass
from scipy.ndimage.interpolation import rotate, shift
import freestream
import frzout

def run_cmd(*args):
    """
    Run and log a subprocess.

    """
    cmd = ' '.join(args)
    logging.info('running command: %s', cmd)

    try:
        proc = subprocess.run(
            cmd.split(), check=True,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            universal_newlines=True
        )
        print(proc.stdout)
        print(proc.stderr)
    except subprocess.CalledProcessError as e:
        logging.error(
            'command failed with status %d:\n%s',
            e.returncode, e.output.strip('\n')
        )
        raise
    else:
        logging.debug(
            'command completed successfully:\n%s',
            proc.stdout.strip('\n')
        )
        return proc

def read_text_file(filename):
    """
    Read a text file into a nested list of bytes objects,
    skipping comment lines (#).

    """
    with open(filename, 'rb') as f:
        return [l.split() for l in f if not l.startswith(b'#')]

class Parser(argparse.ArgumentParser):
    """
    ArgumentParser that parses files with 'key = value' lines.

    """
    def __init__(self, *args, fromfile_prefix_chars='@', **kwargs):
        super().__init__(
            *args, fromfile_prefix_chars=fromfile_prefix_chars, **kwargs
        )

    def convert_arg_line_to_args(self, arg_line):
        # split each line on = and prepend prefix chars to first arg so it is
        # parsed as a long option
        args = [i.strip() for i in arg_line.split('=', maxsplit=1)]
        args[0] = 2*self.prefix_chars[0] + args[0]
        return args


parser = Parser(
    usage=''.join('\n  %(prog)s ' + i for i in [
        '[options] <results_file>',
        'checkpoint <checkpoint_file>',
        '-h | --help',
    ]),
    description='''
Run relativistic heavy-ion collision events.

In the first form, run events according to the given options (below) and write
results to binary file <results_file>.
''',
    formatter_class=argparse.RawDescriptionHelpFormatter
)
parser.add_argument(
    '--logfile', type=os.path.abspath, metavar='PATH',
    help='log file (default: stdout)'
)
parser.add_argument(
    '--loglevel', choices={'debug', 'info', 'warning', 'error', 'critical'},
    default='info',
    help='log level (default: %(default)s)'
)
parser.add_argument(
    'results', type=os.path.abspath,
    help=argparse.SUPPRESS
)
parser.add_argument(
    '--nevents', type=int, metavar='INT',
    help='number of pythia events per pThat bin'
)
parser.add_argument(
    '--nhardevents', type=int, default=20, metavar='INT',
    help='number of pythia events per pThat bin'
)
parser.add_argument(
    '--avg-ic', default='off', metavar='VAR',
    help='if on, generate 500 IC events in the centrality bin and take average'
)
parser.add_argument(
    '--rankvar', metavar='VAR',
    help='environment variable containing process rank'
)
parser.add_argument(
    '--rankfmt', metavar='FMT',
    help='format string for rank integer'
)
parser.add_argument(
    '--tmpdir', type=os.path.abspath, metavar='PATH',
    help='temporary directory (default: {})'.format(tempfile.gettempdir())
)
parser.add_argument(
    '--trento-args', default='', metavar='ARGS',
    help="arguments passed to trento (default: '%(default)s')"
)
parser.add_argument(
    '--lido-args', metavar='ARGS',
    help="arguments passed to lido"
)
parser.add_argument(
    '--tau-fs', type=float, default=.5, metavar='FLOAT',
    help='free streaming time [fm] (default: %(default)s fm)'
)
parser.add_argument(
    '--xi-fs', type=float, default=.5, metavar='FLOAT',
    help='energy loss staring time / tau-fs'
)
parser.add_argument(
    '--hydro-args', default='', metavar='ARGS',
    help='arguments passed to osu-hydro (default: empty)'
)
parser.add_argument(
    '--Tswitch', type=float, default=.150, metavar='FLOAT',
    help='particlization temperature [GeV] (default: %(default).3f GeV)'
)
parser.add_argument(
    '--sqrts', type=float, default=5020, metavar='FLOAT',
    help='center-of-mass energy'
)
parser.add_argument(
    '--proj', type=str, default="Pb", metavar='STR',
    help='projectile'
)
parser.add_argument(
    '--targ', type=str, default="Pb", metavar='STR',
    help='target'
)
parser.add_argument(
    '--setting-path', type=os.path.abspath, metavar='PATH',
    help='lido, pythia, setting path'
)

Light_species= [
            ('pion', 211),
            ('kaon', 321),
            ('proton', 2212),
            ('Lambda', 3122),
            ('Sigma0', 3212),
            ('Xi', 3312),
            ('Omega', 3334)]
Strange_species=[ 
            ('kaon', 321),
            ('Lambda', 3122),
            ('Sigma0', 3212),
            ('Xi', 3312),
            ('Omega', 3334)]
Charm_species = [
            ('D', 411),
            ('D0', 421),
            ('D*+', 413),
            ('D*0', 423)]
Bottom_species=[
            ('B+', 511),
            ('B0', 521),
            ('B*+', 513),
            ('B*0', 523)]

# fully specify numeric data types, including endianness and size, to
# ensure consistency across all machines
float_t = '<f8'
int_t = '<i8'
complex_t = '<c16'
jetRs = np.array([0.2,0.3,0.4,0.6,0.8,1.0])
JetRaapTbins = np.array([4,6,8,12,16,20,25,30,35,
                      40,45,50,60,70,80,90,100,
                      110,120,140,160,180,200,
                      240,280,320,360,400,500,
                      600,800,1000,1200,1400,1600,2000,2500])
JetShapepTbins = np.array([10,20,30,40,60,80,100,120,2500])
JetFragpTbins = np.array([10,20,30,40,60,70,100,126,158,200,251,316,398,600,800])
JetShaperbins = np.array([0., .05, .1, .15,  .2, .25, .3,
                         .35, .4, .45, .5,  .6, .7,  .8,
                          1., 1.5, 2.0, 2.5, 3.0])
JetFragzbins = np.array([0, 0.001, 0.002, 
                         0.0035,.005,.0065,.0085, .01,
                        .012, .016, .02,.025, .03,.04,.05,
                        .07, .09, .120, .15, .20, 
                        .25, .35, .45, .55, .65, .8,
                        1.])
JetFragzpTbins = np.array([
         0.0, 0.15, 0.3, 0.5       ,   0.70626877,   0.99763116,   1.40919147,
         1.99053585,   2.81170663,   3.97164117,   5.61009227,
         7.92446596,  11.19360569,  15.8113883 ,  22.33417961,
         31.54786722,  44.56254691,  62.94627059,  88.9139705 ,
         125.59432158, 177.40669462, 250.59361681, 353.97289219,
         500.])
xJpTbins = np.array([10,20,30,40,60,80,100,126,
                     158,178,200,224,251,282,316,398,562])
xJbins = np.array([0, .05,  .1, .15,  .2,
                 .25,  .3, .35,  .4, .45,
                  .5, .55,  .6, .65,  .7,
                 .75,  .8, .85,  .9, .95, 1.0])
HadronRaapTbins = np.array([0, 1.,2.,3,
                            4,5,6,7,8,10,12,14,16,18,
                            20,22,24,28,32,36,40,45,50,
                            55,60,65,70,80,90,100,
                            110,120,140,160,180,200,
                            220,240,260,300,350,
                            400,500,600,800,1000]) 
# results "array" (one element)
result_dtype=[
    ('initial_entropy', float_t),
    ('Ncoll', float_t),
    ('dNchdeta', float_t),
    ('dETdeta', float_t),
    ('Soft-Qn',   [('M', int_t), ('Qn', complex_t, 8)]),
    *[  (s, [('spectra', float_t, [jetRs.shape[0], JetRaapTbins.shape[0]-1]),
             ('vn', [ ('M', float_t, [jetRs.shape[0], JetRaapTbins.shape[0]-1]), 
                      ('Qn', complex_t, [jetRs.shape[0], 2, JetRaapTbins.shape[0]-1])
                    ]),
             ('shape', [ ('dsigma', float_t, JetShapepTbins.shape[0]-1),
                         ('dpTdr', float_t, [JetShapepTbins.shape[0]-1, 
                                              JetShaperbins.shape[0]-1]) ]
             ),
             ('frag', [ ('dsigma', float_t, JetFragpTbins.shape[0]-1),
                        ('dNdz', float_t, [JetFragpTbins.shape[0]-1, 
                                           JetFragzbins.shape[0]-1]),
                        ('dNdpT', float_t, [JetFragpTbins.shape[0]-1, 
                                            JetFragzpTbins.shape[0]-1]) ])
           ]) for s in ['Inclusive-jets', 'c-jets', 'b-jets']
     ],
    *[ (s, [  ('shape', [ ('dsigma-low-pT', float_t, 1),
                        ('dsigma-high-pT', float_t, 1),
                      ('dNdr-low-pT', float_t, JetShaperbins.shape[0]-1),
                      ('dNdr-high-pT', float_t, JetShaperbins.shape[0]-1) ]),
              ('frag',  [ ('dsigma', float_t,  JetFragpTbins.shape[0]-1),
                          ('dNdz', float_t, [JetFragpTbins.shape[0]-1,
                                             JetFragzbins.shape[0]-1]) ])
           ]) for s in ['c-jet-corr', 'b-jet-corr']
    ],
    *[ (s, [ ('spectra', float_t, HadronRaapTbins.shape[0]-1),
             ('vn', [ ('M', float_t, HadronRaapTbins.shape[0]-1),
                      ('Qn', complex_t, [2, HadronRaapTbins.shape[0]-1])
                    ])  
           ]) for s in ['chg','s','c','b']
     ],
    ('xJ', [ ('dsigma', float_t, xJpTbins.shape[0]-1),
             ('dNdxJ', float_t, [xJpTbins.shape[0]-1,xJbins.shape[0]-1])
           ]
    )
]

class StopEvent(Exception):
    """ Raise to end an event early. """

# Run IC
def _initial_conditions(nevents=1, initial_file='initial.hdf', avg='off'):
    def average_ic(fname):
        with h5py.File(fname, 'a') as f:
            densityavg = np.zeros_like(f['event_0/matter_density'].value)
            Ncollavg = np.zeros_like(f['event_0/Ncoll_density'].value)
            dxy = f['event_0'].attrs['dxy']
            Neve = len(f.values())
            for eve in f.values():
                # step1, center the event
                NL = int(eve.attrs['Nx']/2)
                density = eve['matter_density'].value
                comxy = -np.array(center_of_mass(density))+np.array([NL, NL])
                density = shift(density, comxy)
                Ncoll = shift(eve['Ncoll_density'].value, comxy)
                # step2, rotate the event to align psi2
                psi2 = eve.attrs['psi2']
                imag_psi2 = psi2*180./np.pi + (90. if psi2<0 else -90.)
                densityavg += rotate(density, angle=imag_psi2, reshape=False)
                Ncollavg += rotate(Ncoll, angle=imag_psi2, reshape=False)    
            # step3 take average
            densityavg /= Neve
            Ncollavg /= Neve
        # rewrite the initial.hdf file with average ic
        with h5py.File(fname, 'w') as f:
            gp = f.create_group('avg_event')
            gp.create_dataset('matter_density', data=densityavg)
            gp.create_dataset('Ncoll_density', data=Ncollavg) 
            gp.attrs.create('Nx', densityavg.shape[1])
            gp.attrs.create('Ny', densityavg.shape[0])
            gp.attrs.create('dxy', dxy)            
    try:
        os.remove(initial_file)
    except FileNotFoundError:
        pass    
    if avg == 'on':
        logging.info("averaged initial condition mode, could take a while")
    trentoargs = args.trento_args
    run_cmd(
        'trento',
        '{} {}'.format(args.proj, args.targ),
        '--number-events {}'.format(nevents if avg=='off' else 500),
        '--grid-step {} --grid-max {}'.format(grid_step, grid_max_target),
        '--output', initial_file,
        '{}'.format(trentoargs)
    )
    if avg == 'on':
        logging.info("taking average over 500 trento events")
        average_ic(initial_file)

    ### create iterable initial conditon generator
    with h5py.File(initial_file, 'r') as f:
        for dset in f.values():
            ic = np.array(dset['matter_density'])
            nb = np.array(dset['Ncoll_density'])
            yield np.array([ic, nb])

def save_fs_with_hydro(ic):
    tau_eloss = args.xi_fs*args.tau_fs
    # roll ic by index 1 to match hydro
    #ic = np.roll(np.roll(ic, shift=-1, axis=0), shift=-1, axis=1)
    # use same grid settings as hydro output
    with h5py.File('JetData.h5','a') as f:
        taufs = f['Event'].attrs['Tau0'][0]
        dtau = f['Event'].attrs['dTau'][0]
        dxy = f['Event'].attrs['DX'][0]
        ls = f['Event'].attrs['XH'][0]
        n = 2*ls + 1
        coarse = int(dxy/grid_step+.5)
        # [tau0, tau0+dtau, tau0+2*dtau, ..., taufs - dtau] + hydro steps...
        nsteps = int((taufs-tau_eloss)/dtau)
        tau0 = taufs-dtau*nsteps
        if tau0 < 1e-2: # if tau0 too small, skip the first step
            tau0 += dtau
            nsteps -= 1
        logging.info("taueloss={:1.2f}, tau0={:1.2f}".format(tau_eloss, tau0))
        taus = np.linspace(tau0, taufs-dtau, nsteps)
        # First, rename hydro frames and leave the first few name slots to FS
        event_gp = f['Event']
        for i in range(len(event_gp.keys()))[::-1]:
            old_name = 'Frame_{:04d}'.format(i)
            new_name = 'Frame_{:04d}'.format(i+nsteps)
            event_gp.move(old_name, new_name)
        # Second, overwrite tau0 with FS starting time, and save taufs where
        # FS and hydro is separated
        event_gp.attrs.create('Tau0', [tau0])
        event_gp.attrs.create('TauFS', [taufs])
        # Thrid, fill the first few steps with Freestreaming results
        for itau, tau in enumerate(taus):
            frame = event_gp.create_group('Frame_{:04d}'.format(itau))
            fs = freestream.FreeStreamer(ic, grid_max, tau)
            for fmt, data, arglist in [
                ('e', fs.energy_density, [()]),
                ('V{}', fs.flow_velocity, [(1,), (2,)]),
                ('Pi{}{}', fs.shear_tensor, [(0,0), (0,1), (0,2),
                                                        (1,1), (1,2),
                                                               (2,2)] ),
                ]:
                for a in arglist:
                    X = data(*a).T # to get the correct x-y with vishnew
                    if fmt == 'V{}': # Convert u1, u2 to v1, v2
                        X = X/data(0).T
                    X = X[::coarse, ::coarse]
                    diff = X.shape[0] - n
                    start = int(abs(diff)/2)
                    if diff > 0:
                        # original grid is larger -> cut out middle square
                        s = slice(start, start + n)
                        X = X[s, s]
                    elif diff < 0:
                        # original grid is smaller
                        #  -> create new array and place original grid in middle
                        Xn = np.zeros((n, n))
                        s = slice(start, start + X.shape[0])
                        Xn[s, s] = X
                        X = Xn
                    if fmt == 'V{}':
                        Comp = {1:'x', 2:'y'}
                        frame.create_dataset(fmt.format(Comp[a[0]]), data=X)
                    if fmt == 'e':
                        frame.create_dataset(fmt.format(*a), data=X)
                        frame.create_dataset('P', data=X/3.)
                        frame.create_dataset('BulkPi', data=X*0.)
                        prefactor = 1.0/15.62687/5.068**3 
                        frame.create_dataset('Temp', data=(X*prefactor)**0.25)
                        s = (X + frame['P'].value)/(frame['Temp'].value+1e-14)
                        frame.create_dataset('s', data=s)
                    if fmt == 'Pi{}{}': 
                        frame.create_dataset(fmt.format(*a), data=X)
            pi33 = -(frame['Pi00'].value + frame['Pi11'].value \
                                             + frame['Pi22'].value)
            frame.create_dataset('Pi33', data=pi33)
            pi3Z = np.zeros_like(pi33)
            frame.create_dataset('Pi03', data=pi3Z)
            frame.create_dataset('Pi13', data=pi3Z)
            frame.create_dataset('Pi23', data=pi3Z)

def run_hydro(ic, event_size, coarse=False, dt_ratio=.25):
    # append switching energy density to hydro arguments
    eswitch = hrg.energy_density()
    hydro_args = [args.hydro_args, 'edec={}'.format(eswitch)]
    hydro_args_coarse = [
        'etas_hrg=0 etas_min=0 etas_slope=0 zetas_max=0 zetas_width=0',
        'edec={}'.format(frzout.HRG(.110, **hrg_kwargs).energy_density())
    ]
    # first freestream
    fs = freestream.FreeStreamer(ic, grid_max, args.tau_fs)
    dxy = grid_step * (coarse or 1)
    ls = math.ceil(event_size/dxy)  # the osu-hydro "ls" parameter
    n = 2*ls + 1  # actual number of grid cells
    for fmt, f, arglist in [
            ('ed', fs.energy_density, [()]),
            ('u{}', fs.flow_velocity, [(1,), (2,)]),
            ('pi{}{}', fs.shear_tensor, [(1, 1), (1, 2), (2, 2)]),
    ]:
        for a in arglist:
            X = f(*a)
            if coarse:
                X = X[::coarse, ::coarse]
            diff = X.shape[0] - n
            start = int(abs(diff)/2)
            if diff > 0:
                # original grid is larger -> cut out middle square
                s = slice(start, start + n)
                X = X[s, s]
            elif diff < 0:
                # original grid is smaller
                #  -> create new array and place original grid in middle
                Xn = np.zeros((n, n))
                s = slice(start, start + X.shape[0])
                Xn[s, s] = X
                X = Xn
            X.tofile(fmt.format(*a) + '.dat')
    dt = dxy*dt_ratio
    run_cmd(
        'osu-hydro',
        't0={} dt={} dxy={} nls={}'.format(args.tau_fs, dt, dxy, ls),
        *(hydro_args_coarse if coarse else hydro_args)
    )
    surface = np.fromfile('surface.dat', dtype='f8').reshape(-1, 26)
    # surface columns:
    #   0     1  2  3    
    #   tau  x  y  eta  
    #   4     5     6     7
    #   dsigma_t  dsigma_x  dsigma_y  dsigma_z
    #   8    9    10
    #   v_x  v_y  v_z
    #   11    12    13    14    
    #   pitt  pitx  pity  pitz
    #     15    16    17
    #     pixx  pixy  pixz
    #           18    19
    #           piyy  piyz
    #             20
    #             pizz
    #   21   22   23   24   25
    #   Pi   T    e    P    muB
    if not coarse:
        logging.info("Save free streaming history with hydro histroy")
        save_fs_with_hydro(ic)
    # end event if the surface is empty -- this occurs in ultra-peripheral
    # events where the initial condition doesn't exceed Tswitch
    if surface.size == 0:
        raise StopEvent('empty surface')
    # pack surface data into a dict suitable for passing to frzout.Surface
    return dict(
            x=surface[:, 0:3],
            sigma=surface[:, 4:7],
            v=surface[:, 8:10],
            pi=dict(xx=surface.T[15],xy=surface.T[16], yy=surface.T[18]),
            Pi=surface.T[21]
        )
def save(results):
    # Load jet results and put them into the results file
    # 1 leading particle spectra and vn
    pT, pTl, pTh, dnch, dnpi, dnk, dnD, dnB = \
        np.loadtxt("./lido-LeadingHadron.dat").T
    results['chg']['spectra'] = dnch
    results['s']['spectra'] = dnk
    results['c']['spectra'] = dnD
    results['b']['spectra'] = dnB

    data = np.loadtxt("./lido-LeadingQn.dat").T
    results['chg']['vn']['M'] = data[3]
    results['chg']['vn']['Qn'][0] = data[4]+1j*data[5]
    results['chg']['vn']['Qn'][1] = data[6]+1j*data[7]
    results['s']['vn']['M'] = data[8]
    results['s']['vn']['Qn'][0] = data[9]+1j*data[10]
    results['s']['vn']['Qn'][1] = data[11]+1j*data[12]
    results['c']['vn']['M'] = data[13]
    results['c']['vn']['Qn'][0] = data[14]+1j*data[15]
    results['c']['vn']['Qn'][1] = data[16]+1j*data[17]
    results['b']['vn']['M'] = data[18]
    results['b']['vn']['Qn'][0] = data[19]+1j*data[20]
    results['b']['vn']['Qn'][1] = data[21]+1j*data[22]

    # 2 Jet spectra and vn
    for iR, R in enumerate(jetRs):
        pT, pTl, pTh, dJ, dD, dB = \
            np.loadtxt("./lido-jet-R-{}-spectra.dat".format(iR)).T
        results['Inclusive-jets']['spectra'][iR] = dJ
        results['c-jets']['spectra'][iR] = dD
        results['b-jets']['spectra'][iR] = dB

        data = \
            np.loadtxt("./lido-jet-R-{}-qn.dat".format(iR)).T
        results['Inclusive-jets']['vn']['M'][iR] = data[3]
        results['Inclusive-jets']['vn']['Qn'][iR][0] = data[4]+1j*data[5]
        results['Inclusive-jets']['vn']['Qn'][iR][1] = data[6]+1j*data[7]

    # 3 Jet shape and fragmentation function
    data = np.loadtxt("./lido-jetshape.dat")
    for i, d in enumerate(data):
        results['Inclusive-jets']['shape']['dsigma'] = d[3]
        results['Inclusive-jets']['shape']['dpTdr'] = d[4:]
    data = np.loadtxt("./lido-D-jetshape.dat")
    for i, d in enumerate(data):
        results['c-jets']['shape']['dsigma'] = d[3]
        results['c-jets']['shape']['dpTdr'] = d[4:]
    data = np.loadtxt("./lido-B-jetshape.dat")
    for i, d in enumerate(data):
        results['b-jets']['shape']['dsigma'] = d[3]
        results['b-jets']['shape']['dpTdr'] = d[4:]

    data = np.loadtxt("./lido-jet-dNdz.dat")
    datapT = np.loadtxt("./lido-jet-dNdpT.dat")
    for i, (d, dpT) in enumerate(zip(data,datapT)):
        results['Inclusive-jets']['frag']['dsigma'] = d[3]
        results['Inclusive-jets']['frag']['dNdz'] = d[4:]
        results['Inclusive-jets']['frag']['dNdpT'] = dpT[4:]
    data = np.loadtxt("./lido-D-jet-dNdz.dat")
    datapT = np.loadtxt("./lido-D-jet-dNdpT.dat")
    for i, (d, dpT) in enumerate(zip(data,datapT)):
        results['c-jets']['frag']['dsigma'] = d[3]
        results['c-jets']['frag']['dNdz'] = d[4:]
    data = np.loadtxt("./lido-B-jet-dNdz.dat")
    datapT = np.loadtxt("./lido-B-jet-dNdpT.dat")
    for i, (d, dpT) in enumerate(zip(data,datapT)):
        results['b-jets']['frag']['dsigma'] = d[3]
        results['b-jets']['frag']['dNdz'] = d[4:]   
   
    # 4 Hard-hard correlations: jet-flavor correlation
    data = np.loadtxt("./lido-D-in-jet.dat")
    for i, d in enumerate(data):
        results['c-jet-corr']['frag']['dsigma'] = d[3]
        results['c-jet-corr']['frag']['dNdz'] = d[4:]
    data = np.loadtxt("./lido-B-in-jet.dat")
    for i, d in enumerate(data):
        results['b-jet-corr']['frag']['dsigma'] = d[3]
        results['b-jet-corr']['frag']['dNdz'] = d[4:]

    d = np.loadtxt("./lido-dDdr.dat")
    results['c-jet-corr']['shape']['dsigma-low-pT'] = d[0,3]
    results['c-jet-corr']['shape']['dNdr-low-pT'] = d[0,4:]
    results['c-jet-corr']['shape']['dsigma-high-pT'] = d[1,3]
    results['c-jet-corr']['shape']['dNdr-high-pT'] = d[1,4:]

    d = np.loadtxt("./lido-dBdr.dat")
    results['b-jet-corr']['shape']['dsigma-low-pT'] = d[0,3]
    results['b-jet-corr']['shape']['dNdr-low-pT'] = d[0,4:]
    results['b-jet-corr']['shape']['dsigma-high-pT'] = d[1,3]
    results['b-jet-corr']['shape']['dNdr-high-pT'] = d[1,4:]   
   
    # 5 Hard-hard correlations: di-jet asymmetry
    data = np.loadtxt("./lido-xJ.dat")
    for i, d in enumerate(data):
        results['xJ']['dsigma'] = d[3]
        results['xJ']['dNdxJ'] = d[4:]

def run_pp_event(results, ic, nb, event_number):
    run_cmd(
        'Lido_pp',
        '-y {:s}/pythia_pp.txt'.format(args.setting_path),
        '--ic ./initial.hdf',
        '--eid {:d}'.format(0 if args.nevents is None else event_number-1),
        '--pTtrack .7',
        '-n {:d}'.format(args.nhardevents),
        '--jet',
        '-o ./'
    )
    save(results)

def run_single_event(results, ic, nb, event_number):
    logging.info(
        'free streaming initial condition for %.3f fm',
        args.tau_fs
    )
    fs = freestream.FreeStreamer(ic, grid_max, args.tau_fs)
    # run coarse event on large grid and determine max radius
    rmax = math.sqrt((
        run_hydro(ic, event_size=27, coarse=3)['x'][:, 1:3]**2
        ).sum(axis=1).max())
    logging.info('rmax = %.3f fm', rmax)
    # now run normal event with size set to the max radius
    # and create sampler surface object
    surface = frzout.Surface(**run_hydro(ic, event_size=rmax), ymax=3.5)
    logging.info('%d freeze-out cells', len(surface))
    logging.info('sampling surface with frzout')
    E, px, py, pz = frzout.sample(surface, hrg)['p'].T
    nparts = len(E)
    ET = np.sqrt(E**2-pz**2)
    pT = np.sqrt(px**2+py**2)
    pabs = np.sqrt(pT**2+pz**2)
    y = .5*np.log((E+pz)/(E-pz))
    eta = .5*np.log((pabs+pz)/(pabs-pz))
    phi = np.arctan2(py, px)
    logging.info('produced %d particles from QGP', nparts)
    cut = np.abs(eta)<.5
    results['dNchdeta'] = cut.sum()
    results['dETdeta'] = ET[cut].sum()
    if nparts == 0:
        raise StopEvent('no particles produced')
    # oversample to get a high percision event plane at freezeout
    nloop=0
    ncount = 0
    while ncount < 10**6 and nloop < 10000:
        nloop += 1
        E, px, py, pz = frzout.sample(surface, hrg)['p'].T
        pT = np.sqrt(px**2+py**2)
        pabs = np.sqrt(pT**2+pz**2)
        eta = .5*np.log((pabs+pz)/(pabs-pz))
        phi = np.arctan2(py, px)
        cut = (np.abs(eta) < 2.) & (0.3 < pT)
        ncount += cut.sum()
        results['Soft-Qn']['M'] += cut.sum()
        for n in range(8): # Q1 to Q8
            results['Soft-Qn']['Qn'][n] += np.exp(1j*(n+1)*phi[cut]).sum()

    # Run jet
    run_cmd(
        'Lido2DHydro',
        '-y {:s}/pythia.txt'.format(args.setting_path),
        '--ic ./initial.hdf',
        '--eid {:d}'.format(0 if args.nevents is None else event_number-1),
        '--hydro ./JetData.h5',
        '-s {:s}/lido_settings.xml'.format(args.setting_path),
        '-t {:s}/table.h5'.format(args.setting_path),
        '-r {:s}/response.h5'.format(args.setting_path),
        '-n {:d}'.format(args.nhardevents),
        args.lido_args,
        '-o ./'
    )

    save(results)

def run_events(results_file):
    ### 1. Generate initial conditain chain
    # if nevents was specified, generate that number of initial conditions
    # otherwise generate indefinitely
    initial_conditions = (
            chain.from_iterable(_initial_conditions() \
            for _ in repeat(None))
            if args.nevents is None else
            _initial_conditions(args.nevents, avg=args.avg_ic)
        )
    ### 2: event recorder
    results = np.empty((), dtype=result_dtype)    
    nfail = 0

    ### 3: run each initial condition event and save results to file
    for n, (ic, nb) in enumerate(initial_conditions, start=1):
        results.fill(0)
        results['initial_entropy'] = ic.sum() * grid_step**2
        results['Ncoll'] = nb.sum() * grid_step**2
        logging.info("Nb %d", results['Ncoll'])
        logging.info('starting event %d', n)
        try:
            if args.proj == 'p' and args.targ == 'p':
                run_pp_event(results, ic, nb, n)
            else:
                run_single_event(results, ic, nb, n)
        except StopEvent as e:
            logging.info('event stopped: %s', e)
        except Exception:
            logging.exception('event %d failed', n)
            nfail += 1
            if nfail > 3 and nfail/n > .5:
                logging.critical('too many failures, stopping events')
                break
            logging.warning('continuing to next event')
            continue

        results_file.write(results.tobytes())
        logging.info('event %d completed successfully', n)
    return n > nfail

def parse_args():
    def usage():
        parser.print_usage(sys.stderr)
        sys.exit(2)
    if len(sys.argv) == 1:
        usage()
    return parser.parse_args()


def main():
    global grid_step, grid_max_target, grid_n, grid_max, args, hrg_kwargs, hrg
    grid_step = .15
    grid_max_target = 15
    grid_n = math.ceil(2*grid_max_target/grid_step)
    grid_max = .5*grid_n*grid_step
    logging.info(
        'grid step = %.6f fm, n = %d, max = %.6f fm',
        grid_step, grid_n, grid_max
    )
    args = parse_args()
    hrg_kwargs = dict(species='urqmd', res_width=True)
    hrg = frzout.HRG(args.Tswitch, **hrg_kwargs)
    filemode = 'w'
    # must handle rank first since it affects paths
    if args.rankvar:
        rank = os.getenv(args.rankvar)
        if rank is None:
            sys.exit('rank variable {} is not set'.format(args.rankvar))
        if args.rankfmt:
            rank = args.rankfmt.format(int(rank))

        # append rank to path arguments, e.g.:
        #   /path/to/output.log  ->  /path/to/output/<rank>.log
        for a in ['results']:
            value = getattr(args, a)
            if value is not None:
                root, ext = os.path.splitext(value)
                setattr(args, a, os.path.join(root, rank) + ext)

    os.makedirs(os.path.dirname(args.results), exist_ok=True)
    if args.logfile is None:
        logfile_kwargs = dict(stream=sys.stdout)
    else:
        logfile_kwargs = dict(filename=args.logfile, filemode=filemode)
        os.makedirs(os.path.dirname(args.logfile), exist_ok=True)

    logging.basicConfig(
        level=getattr(logging, args.loglevel.upper()),
        format='[%(levelname)s@%(relativeCreated)d] %(message)s',
        **logfile_kwargs
    )
    logging.captureWarnings(True)

    start = datetime.datetime.now()
    logging.info('started at %s', start)
    logging.info('arguments: %r', args)

    # translate SIGTERM to KeyboardInterrupt
    signal.signal(signal.SIGTERM, signal.default_int_handler)
    logging.debug('set SIGTERM handler')

    with open(args.results, filemode + 'b') as results_file, \
         tempfile.TemporaryDirectory(prefix='hic-', dir=args.tmpdir) as workdir:
        os.chdir(workdir)
        logging.info('working directory: %s', workdir)
        try:
            status = run_events(results_file)
        except KeyboardInterrupt:
            # after catching the initial SIGTERM or interrupt, ignore them
            # during shutdown -- this ensures everything will exit gracefully
            # in case of additional signals (short of SIGKILL)
            signal.signal(signal.SIGTERM, signal.SIG_IGN)
            signal.signal(signal.SIGINT, signal.SIG_IGN)
            status = True
            logging.info(
                'interrupt or signal at %s, cleaning up...',
                datetime.datetime.now()
            )
    end = datetime.datetime.now()
    logging.info('finished at %s, %s elapsed', end, end - start)

    if not status:
        sys.exit(1)

if __name__ == "__main__":
    main()

