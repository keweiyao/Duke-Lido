import subprocess
import sys
import os
from multiprocessing import Pool, cpu_count
from subprocess import Popen, PIPE

# def call_JS(pTmin, pTmax):
#    subprocess.call('./LidoJetTest {:d} {:d} {:d}'.format(pTmin, pTmax, 1000), shell=True)
#Nproc = cpu_count()
#print("CPU COUNT: "+ str(Nproc))
#pTbins = [{30, 50, 70, 90, 110, 130, 150, 200, 300, 400, 600}]
#pool = Pool(Nproc)
#pool.starmap(call_JS, zip(pTbins[:-1], pTbins[1:]))


##################################################################################

os.makedirs("./results", exist_ok=True)

pTbins = [10, 20, 30, 50, 70, 90, 110, 150, 200, 300, 500, 700]
# for i in range(len(pTbins)-1):
#    subprocess.call("cp jetscape_init.xml ./group_input/jetscape_init_{:d}.xml".format(i), shell=True)
# jet-pythia -y ../examples/pythia-jet-setting -s settings.xml -t ../table.h5 --pthat-low 150 --pthat-high 200

subprocess.run("../build/Lido-TabGen -s ../examples/lido_settings.xml -t ../build/table.h5", shell=True, check=True)
#process = Popen(['../build/Lido-TabGen', "-s", "../examples/settings.xml",
#                 "-t", "../build/table.h5"], stdout=PIPE, stderr=PIPE)
#stdout, stderr = process.communicate()
#process.wait()

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
           "-n", "10000"]
    return arg


procs = [Popen([cmd,
                *getarg(pTbins[i], pTbins[i+1])
                ]
               )
         for i in range(len(pTbins)-1)]
for p in procs:
    p.wait()
