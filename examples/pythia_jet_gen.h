#ifndef PYTHIA_WRAPPER_H
#define PYTHIA_WRAPPER_H

#include "simpleLogger.h"
#include "Pythia8/Pythia.h"
#include "workflow.h"
#include <sstream>

using namespace Pythia8;

class PythiaGen{
public: 
    PythiaGen(std::string f_pythia, int pTHL, int pTHH);
    void Generate(std::vector<particle> & plist);
    double sigma_gen(void){
        return pythia.info.sigmaGen();
    }
private:
    Pythia pythia;
};

PythiaGen::PythiaGen(std::string f_pythia, int pTHL, int pTHH) {    
    // read pythia settings
    pythia.readFile(f_pythia);
    // suppress output
    pythia.readString("SoftQCD:all = off");
    pythia.readString("PromptPhoton:all=off");
    pythia.readString("WeakSingleBoson:all=off");
    pythia.readString("WeakDoubleBoson:all=off");
    pythia.readString("Init:showProcesses = off");  
    pythia.readString("Init:showMultipartonInteractions = off");  
    pythia.readString("Init:showChangedSettings = off");  
    pythia.readString("Init:showChangedParticleData = off");  
    pythia.readString("Next:numberCount = 1000");  
    pythia.readString("Next:numberShowInfo = 0");  
    pythia.readString("Next:numberShowProcess = 0");  
    pythia.readString("Next:numberShowEvent = 0"); 

    std::ostringstream s1, s2;
    s1 << "PhaseSpace:pTHatMin = " << pTHL;
    s2 << "PhaseSpace:pTHatMax = " << pTHH;
    pythia.readString(s1.str());
    pythia.readString(s2.str());
    // Init
    pythia.init();
}

void PythiaGen::Generate(std::vector<particle> & plist){
    plist.clear();
    pythia.next();
    double weight = pythia.info.sigmaGen();
    for (size_t i = 0; i < pythia.event.size(); ++i) 
    {
        auto p = pythia.event[i];
        if (p.isFinal()) {
            // final momenta 
            fourvec p0{p.e(), p.px(), p.py(), p.pz()};
            particle _p;
            _p.pid = p.id(); 
            _p.mass = std::abs(p.m());
	        _p.x0 = fourvec{0,0,0,0}; 
	        _p.x = _p.x0; 
	        _p.p0 = p0; 
	        _p.p = _p.p0; 
            _p.weight = weight;
            _p.is_vac = false;
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
