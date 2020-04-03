#ifndef PYTHIA_WRAPPER_H
#define PYTHIA_WRAPPER_H

#include "simpleLogger.h"
#include "Pythia8/Pythia.h"
#include "workflow.h"
#include <sstream>
#include "predefine.h"

using namespace Pythia8;

class PythiaGen{
public: 
    PythiaGen(std::string f_pythia, std::string f_trento, 
              int pTHL, int pTHH, int iev, double _Q0);
    void Generate(std::vector<particle> & plist);
    double sigma_gen(void){
        return sigma0;
    }
private:
    Pythia pythia, pythia2;
    TransverPositionSampler TRENToSampler;
    double sigma0, Q0;
};

PythiaGen::PythiaGen(std::string f_pythia, std::string f_trento,
                     int pTHL, int pTHH, int iev, double _Q0):
TRENToSampler(f_trento, iev)
{    
    Q0 = _Q0;
    // read pythia settings
    pythia.readFile(f_pythia);
    // suppress output
    pythia.readString("Print:quiet = off");
    pythia.readString("SoftQCD:all = off");
    pythia.readString("PromptPhoton:all=off");
    pythia.readString("WeakSingleBoson:all=off");
    pythia.readString("WeakDoubleBoson:all=off");
    pythia.readString("SpaceShower:QEDshowerByQ=off");
    pythia.readString("TimeShower:QEDshowerByQ = off");
    pythia.readString("Init:showProcesses = off");  
    pythia.readString("Init:showMultipartonInteractions = off");  
    pythia.readString("Init:showChangedSettings = off");  
    pythia.readString("Init:showChangedParticleData = off");  
    pythia.readString("Next:numberCount = 1000");  
    pythia.readString("Next:numberShowInfo = 0");  
    pythia.readString("Next:numberShowProcess = 0");  
    pythia.readString("Next:numberShowEvent = 0"); 

    int processid = getpid();

    std::ostringstream s1, s2, s3, s4;
    s1 << "PhaseSpace:pTHatMin = " << pTHL;
    s2 << "PhaseSpace:pTHatMax = " << pTHH;
    s3 << "Random:seed = " << processid;
    s4 << "TimeShower:pTmin = " << Q0;
    pythia.readString(s1.str());
    pythia.readString(s2.str());
    pythia.readString(s3.str());
    pythia.readString(s4.str());
    // Init
    pythia.init();
    for (int i=0; i<1000; i++) pythia.next();
    sigma0 = pythia.info.sigmaGen();
}

void PythiaGen::Generate(std::vector<particle> & plist){
    double x, y;
    TRENToSampler.SampleXY(x, y);
    plist.clear();
    pythia.next();
    auto & event = pythia.event;
    color_count = event.lastColTag()+1;
    for (size_t i = 0; i < event.size(); ++i) {
        auto p = event[i];
        if (p.isFinal()) {
            LOG_INFO << p.id() << " " << p.col() << " " << p.acol() << " " << color_count;
            // final momenta 
            fourvec p0{p.e(), p.px(), p.py(), p.pz()};
            particle _p; 
            _p.pid = p.id();
            _p.mass = std::abs(p.m());
            _p.x0 = fourvec{0,x,y,0};
            _p.x = _p.x0; 
            _p.tau_i = 0.;
            _p.p0 = p0;
            _p.Q0 = Q0;
 
            if (std::abs(_p.pid) != 4 && 
                        std::abs(_p.pid) != 5 && p.isParton()) {
                _p.mass = 0;
            }
            _p.col = p.col();
            _p.acol = p.acol();
            _p.p0.a[0] = std::sqrt(_p.p0.pabs2()+_p.mass*_p.mass);
            _p.p = _p.p0; 
            _p.weight = sigma0;
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
}


#endif
