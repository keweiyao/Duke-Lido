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

os.makedirs("./results", exist_ok=True)

pTbins=[10, 20, 30, 50, 70, 90, 110, 150, 200, 300, 500, 700]
#for i in range(len(pTbins)-1):
#    subprocess.call("cp jetscape_init.xml ./group_input/jetscape_init_{:d}.xml".format(i), shell=True)
# jet-pythia -y ../examples/pythia-jet-setting -s settings.xml -t ../table.h5 --pthat-low 150 --pthat-high 200

cmd = "../build/jet-pythia" 
def getarg(ptl, pth):
    arg = ["-y", "../examples/pythia-jet-setting",
       "-s", "../examples/settings.xml", 
       "-t", "../table.h5",
       "--tf", "5.0",
       "--temp", "0.35",
       "--afix", "0.3",
       "--pthat-low", "{:d}".format(ptl),
       "--pthat-high", "{:d}".format(pth),
       "-n", "1000"]
    return arg
procs = [ Popen( [ cmd, 
                   *getarg(pTbins[i], pTbins[i+1]) 
                 ] \
                ) \
         for i in range(len(pTbins)-1) ]
for p in procs:
    p.wait()
