#ifndef PYTHIA_WRAPPER_H
#define PYTHIA_WRAPPER_H

#include "simpleLogger.h"
#include "Pythia8/Pythia.h"
#include "workflow.h"
#include <sstream>
#include "predefine.h"
#include "random.h"

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

void reference_pmu(int i, double & tau, Event & event){
    auto p = event[i];
    int im1 = p.mother1();
    int im2 = p.mother2();
    if (im1==im2 && im2 > 0){
        // simply recoil effect, go on
        reference_pmu(im1, tau, event);
    }
    if (im1==0 && im2 ==0) {
        tau = 0.;
	return;
    }
    if (im1>0 && im2==0){
	// radiation
	auto P = event[im1];
	auto k = event[P.daughter1()];
	auto q = event[P.daughter2()];
	fourvec Pmu{P.e(), P.px(), P.py(), P.pz()};
        fourvec kmu{k.e(), k.px(), k.py(), k.pz()};
        fourvec qmu{q.e(), q.px(), q.py(), q.pz()};
	double xk = kmu.t()/(kmu.t()+qmu.t());
	double xq = 1.-xk;
	double Mk2 = k.m()*k.m();
	double Mq2 = q.m()*q.m();
        double MP2 = P.m()*P.m();
	double kT2 = measure_perp(kmu+qmu, kmu).pabs2();
	double tauf = 2*xq*xk*(kmu.t()+qmu.t())/(kT2 + xq*Mk2 + xk*Mq2 - xk*xq*MP2);
	if (p.e() > 0.5*(kmu.t()+qmu.t())){
            tau += 0.;
            reference_pmu(im1, tau, event);
	}
	else {
	    tau += tauf;
	    reference_pmu(im1, tau, event);
	}
    }
    if (im1 != im2 && im1 > 0 && im2 > 0){
        // hard products, reference particles
        return;
    }
}

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
    std::cout<< s4.str();
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
            // final momenta 
            fourvec p0{p.e(), p.px(), p.py(), p.pz()};
            particle _p; 
            _p.pid = p.id();
            _p.mass = std::abs(p.m());
            _p.x0 = fourvec{0,x,y,0};
            _p.x = _p.x0; 
	    double t0 =  0.;
	    reference_pmu(i, t0, event);
            _p.tau_i = t0;
	    LOG_INFO << p0 << " at " << t0/5.076 << " or " << t0/5.076*p0.xT()/p0.t();
            _p.p0 = p0;
            _p.Q0 = Q0;
            _p.Q00 = Q0;
 
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
            _p.charged = p.isCharged(); 
            plist.push_back(_p);
        }
    }
}


#endif
