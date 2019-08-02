import subprocess, sys
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

Rlist=[0.2,0.4,0.6,0.8,1.0]
procs = [ Popen(["./jet_finding", str(Rlist[i]), str(0), str(2.8), str(1000)]) for i in range(len(Rlist)) ]
for p in procs:
   p.wait()
