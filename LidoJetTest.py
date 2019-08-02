import subprocess, sys, os
from multiprocessing import Pool, cpu_count
from subprocess import Popen

#def call_JS(pTmin, pTmax):
#    subprocess.call('./LidoJetTest {:d} {:d} {:d}'.format(pTmin, pTmax, 1000), shell=True)	
#Nproc = cpu_count()
#print("CPU COUNT: "+ str(Nproc))
#pTbins = [{30, 50, 70, 90, 110, 130, 150, 200, 300, 400, 600}]
#pool = Pool(Nproc)
#pool.starmap(call_JS, zip(pTbins[:-1], pTbins[1:]))


##################################################################################

#os.makedirs("group_input", exist_ok=True)

pTbins=[10, 20, 30, 50, 70, 90, 110, 150, 200, 300, 500, 700]
#for i in range(len(pTbins)-1):
#    subprocess.call("cp jetscape_init.xml ./group_input/jetscape_init_{:d}.xml".format(i), shell=True)

procs = [ Popen(["./LidoJetTest", "./group_input/jetscape_init_{:d}.xml".format(i) , str(pTbins[i]), str(pTbins[i+1]), str(1000)]) for i in range(len(pTbins)-1) ]
for p in procs:
    p.wait()
