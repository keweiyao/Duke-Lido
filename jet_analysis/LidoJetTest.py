import subprocess
import sys
import os
from multiprocessing import Pool, cpu_count
from subprocess import Popen, PIPE

##################################################################################

os.makedirs("./results", exist_ok=True)

pTbins = [10, 20, 30, 50, 70, 90, 110, 150, 200, 300, 500, 700]

#generate the table first
subprocess.run("../build/Lido-TabGen -s ../examples/lido_settings.xml -t ../build/table.h5", shell=True, check=True)

cmd = "../build/jet-hydro-couple"

def getarg(ptl, pth):
    arg = ["-y", "../examples/pythia-jet-setting",
           "-s", "../examples/lido_settings.xml",
           "-t", "../build/table.h5",
           "--hydro", "../build/JetData.h5",
           "--ic", "../build/ic.hdf5",
           "--eid", "0",
           "--pthat-low", "{:d}".format(ptl),
           "--pthat-high", "{:d}".format(pth),
           "-n", "2000"]
    return arg

procs = [Popen([cmd, *getarg(pTbins[i], pTbins[i+1])]) for i in range(len(pTbins)-1)]
for p in procs:
    p.wait()
