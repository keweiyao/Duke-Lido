import subprocess
import sys
import os
from multiprocessing import Pool, cpu_count
from subprocess import Popen
from shutil import copyfile

##################################################################################

os.makedirs("./jet_results", exist_ok=True)
#store the current lido_settings and pythia_settings as well
copyfile("../examples/lido_settings.xml", "./jet_results/lido_settings.xml")
copyfile("../examples/pythia-jet-setting.txt", "./jet_results/pythia-jet-setting.txt")

Rlist = [0.2, 0.4, 0.6, 0.8, 1.0]
procs = [Popen(["./jet_finding", str(5020), str(Rlist[i]),
                str(0), str(2.8), str(200)]) for i in range(len(Rlist))]
for p in procs:
    p.wait()
