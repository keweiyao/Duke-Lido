#ifndef PYTHIA_WRAPPER_H
#define PYTHIA_WRAPPER_H

#include "simpleLogger.h"
#include "Pythia8/Pythia.h"
#include "workflow.h"
#include <sstream>
#include "predefine.h"
#include "random.h"


class PGun{
public: 
    PGun(std::string fname);
    void Generate(int pid, int N, double pT, std::vector<particle> & plist);
private:
    TransverPositionSampler TRENToSampler;
};

PGun::PGun(std::string f_trento):TRENToSampler(f_trento, 0){   
    int processid = getpid();
}

void PGun::Generate(int pid, int N, double pT, std::vector<particle> & plist){
    plist.clear();
    double M = 0;
    for (int i=0; i<N; i++){
        double x, y;
        TRENToSampler.SampleXY(y, x);
            double phi = Srandom::dist_phi(Srandom::gen);
            fourvec p0{std::sqrt(pT*pT+M*M), pT*std::cos(phi), pT*std::sin(phi), 0};
            particle _p; 
            _p.pid = pid;
            _p.mass = M;
            _p.x0 = fourvec{0,x,y,0};
            _p.x = _p.x0; 
            _p.tau_i = 0.;
            _p.p0 = p0;
            _p.p = p0;
            _p.Q0 = 1.;
            _p.Q00 = 1.;
            _p.weight = 1;
            _p.is_virtual = false;
            _p.T0 = 0.;
            _p.Tf = 0.;
            _p.mfp0 = 0.;
            _p.vcell.resize(3);
            _p.vcell[0] = 0.; 
            _p.vcell[1] = 0.; 
            _p.vcell[2] = 0.; 
            _p.radlist.clear();
            plist.push_back(_p);
	}
}


#endif
