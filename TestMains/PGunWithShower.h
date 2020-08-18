#ifndef PGUN_W_SHOWER_H
#define PGUN_W_SHOWER_H

#include "simpleLogger.h"
#include "Pythia8/Pythia.h"
#include "workflow.h"
#include <sstream>

using namespace Pythia8;

class PGunWShower{
public: 
    PGunWShower(double _Q0);
    void Generate(double pT, std::vector<particle> & plist);
private:
    Pythia pythia;
    double Q0;
};

PGunWShower::PGunWShower(double _Q0){    
    Q0 = _Q0;
    // suppress output
    pythia.readString("Print:quiet = off");
    pythia.readString("SoftQCD:all = off");
    pythia.readString("PromptPhoton:all=off");
    pythia.readString("WeakSingleBoson:all=off");
    pythia.readString("WeakDoubleBoson:all=off");
    pythia.readString("SpaceShower:QEDshowerByQ=off");
    pythia.readString("TimeShower:QEDshowerByQ = off");
    pythia.readString("ProcessLevel:all = off");
    pythia.readString("HadronLevel:all = off");
    pythia.readString("HadronLevel:Decay = off");
    pythia.readString("Next:numberCount = 1000");  
    pythia.readString("Next:numberShowInfo = 0");  
    pythia.readString("Next:numberShowProcess = 0");  
    pythia.readString("Next:numberShowEvent = 0"); 
    pythia.readString("1:m0 = 0");
    pythia.readString("2:m0 = 0");
    pythia.readString("3:m0 = 0");
    pythia.readString("4:m0 = 1.3");
    pythia.readString("5:m0 = 4.2");

    int processid = getpid();
    std::ostringstream s1, s2;
    s1 << "Random:seed = " << processid;
    s2 << "TimeShower:pTmin = " << Q0;
    pythia.readString(s1.str());
    pythia.readString(s2.str());
    pythia.init();
}



void PGunWShower::Generate(double pT, std::vector<particle> & plist){
    plist.clear();
    for (int repeat=0; repeat<1000; repeat++){
        pythia.event.reset();
        pythia.event.append(21, 23, 101, 102, 
                         pT, 0., 0., pT, 0.); 
        pythia.event.append(21, 23, 102, 101, 
                         -pT, 0., 0., pT, 0.); 
        pythia.event[1].scale(pT);
        pythia.event[2].scale(pT);
        pythia.forceTimeShower(1, 2, pT);
	pythia.next();
        auto & event = pythia.event;
	for (size_t i = 0; i < event.size(); ++i) {
		auto p = event[i];
		if (p.isParton() && event[p.daughter1()].isHadron()) {
		    // final momenta 
		    fourvec p0{p.e(), p.px(), p.py(), p.pz()};
		    particle _p; 
		    _p.pid = p.id();
		    _p.mass = std::abs(p.m());
		    _p.x0 = fourvec{0,0,0,0};
		    _p.x = _p.x0; 
		    _p.tau_i = 0.;
		    _p.p0 = p0;
		    _p.Q0 = Q0;
		    if (std::abs(_p.pid) != 4 && std::abs(_p.pid) != 5)
		        _p.mass = 0;
                    if (std::abs(_p.pid)==4) _p.mass = 1.3;
                    if (std::abs(_p.pid)==5) _p.mass = 4.2;
		    _p.p0.a[0] = std::sqrt(_p.p0.pabs2()+_p.mass*_p.mass);
		    _p.p = _p.p0; 
		    _p.weight = 1.;
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
}


#endif
