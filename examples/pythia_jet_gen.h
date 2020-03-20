#ifndef PYTHIA_WRAPPER_H
#define PYTHIA_WRAPPER_H

#include "simpleLogger.h"
#include "Pythia8/Pythia.h"
#include "workflow.h"
#include <sstream>

using namespace Pythia8;

class PythiaGen{
public: 
    PythiaGen(std::string f_pythia, std::string f_trento, 
              int pTHL, int pTHH, int iev, double Q0);
    void Generate(std::vector<particle> & plist, int heavyid);
    double sigma_gen(void){
        return sigma0;
    }
private:
    Pythia pythia, pythia2;
    TransverPositionSampler TRENToSampler;
    double sigma0;
};

PythiaGen::PythiaGen(std::string f_pythia, std::string f_trento,
                     int pTHL, int pTHH, int iev, double Q0):
TRENToSampler(f_trento, iev)
{    
    // read pythia settings
    pythia.readFile(f_pythia);
    // suppress output
    pythia.readString("Print:quiet = off");
    pythia.readString("SoftQCD:all = off");
    pythia.readString("PromptPhoton:all=off");
    pythia.readString("WeakSingleBoson:all=off");
    pythia.readString("WeakDoubleBoson:all=off");
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
  
    /*pythia2.readString("Random:setSeed = on");
    pythia2.readString("Random:seed = 0");
    pythia2.readString("ProcessLevel:all = off");
    pythia2.readString("TimeShower:QEDshowerByQ = off");
    pythia2.readString("HadronLevel:all = off");
    pythia2.readString("HadronLevel:decay = off");
    pythia2.init();*/
}


void PythiaGen::Generate(std::vector<particle> & plist, int heavyid){
    double x, y;
    TRENToSampler.SampleXY(x, y);
    //LOG_INFO << x/5.076 << " " << y/5.076; 
    bool triggered = false;
    do{
        if (heavyid < 4) triggered = true;
        plist.clear();
        pythia.next();
        ////// Force time shower
        /*pythia2.event.reset();
        // copy old event:
        int k = 1;
        for (int i = 0; i < pythia.event.size(); ++i){
            auto p = pythia.event[i];
            if (p.isFinal() && p.isParton()){
                pythia2.event.append(p.id(), p.status(), p.col(), p.acol(), 
                                     p.px(), p.py(), p.pz(), p.e(), p.m());
                pythia2.event[k].scale(p.scale());
                k++;
            }
        }
        pythia2.forceTimeShower(1, k, .4);
        pythia2.next();*/
        ////////////////////////
        double weight = sigma0;
        auto & event = pythia.event;
	for (size_t i = 0; i < event.size(); ++i) {
		auto p = event[i];
                //if (p.isParton() && event[p.daughter1()].isHadron()
                //    && std::abs(p.y()) < 4) {
		if (p.isFinal() && std::abs(p.y())< 4) {
		    // final momenta 
		    fourvec p0{p.e(), p.px(), p.py(), p.pz()};
		    particle _p; 
		    _p.pid = p.id();
		    _p.mass = std::abs(p.m());
		    _p.x0 = fourvec{0,x,y,0};
                    _p.x = _p.x0; 
                    _p.tau_i = 0.;
		    _p.p0 = p0;
 
		    if (std::abs(_p.pid) != 4 && std::abs(_p.pid) != 5) {
		        _p.mass = 0;
		    }
                    int absid = std::abs(_p.pid);
                    if (heavyid==4){
                        if ( absid == 4 || 
                             absid == 411 || absid == 421 ||
                             absid == 413 || absid == 423  
                           ) triggered = true;
                    }
                    if (heavyid==5){
                        if ( absid == 5 || 
                             absid == 511 || absid == 521 ||
                             absid == 513 || absid == 523  
                           ) triggered = true;
                    }
                    _p.p0.a[0] = std::sqrt(_p.p0.pabs2()+_p.mass*_p.mass);
		    _p.p = _p.p0; 
		    _p.weight = weight;
		    _p.is_vac = false;
		    _p.is_virtual = false;
		    _p.is_recoil = false;
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
    } while(!triggered);
}


#endif
