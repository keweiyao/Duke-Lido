#ifndef PYTHIA_WRAPPER_H
#define PYTHIA_WRAPPER_H

#include "simpleLogger.h"
#include "Pythia8/Pythia.h"
#include "workflow.h"
#include <sstream>
#include "predefine.h"
#include "random.h"

using namespace Pythia8;

class PGunWShower{
public: 
    PGunWShower(double _Q0, std::string fname);
    void Generate(double pT, std::vector<particle> & plist);
private:
    TransverPositionSampler TRENToSampler;
    Pythia pythia;
    double Q0;
};

void reference_pmu(int i, double & tau, Event & event){
    auto p = event[i];
    int absid = p.idAbs();
    int im1 = p.mother1();
    int im2 = p.mother2();
    if (im1==im2 && im2 > 0){
        // simply recoil effect, go on
        reference_pmu(im1, tau, event);
    }
    if (im1==0 && im2 ==0) {
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
            reference_pmu(im1, tau, event);
	}
	else {
	    tau += tauf;
            reference_pmu(im1, tau, event);
	}      
    }
    if (im1 != im2 && im1 > 0 && im2 > 0){
        return;
    }
}

PGunWShower::PGunWShower(double _Q0, std::string f_trento):TRENToSampler(f_trento, 0){   
    Q0 = _Q0;
    pythia.readString("Print:quiet = off");
    pythia.readString("SoftQCD:all = off");
    pythia.readString("PromptPhoton:all=off");
    pythia.readString("WeakSingleBoson:all=off");
    pythia.readString("WeakDoubleBoson:all=off");
    pythia.readString("SpaceShower:QEDshowerByQ=off");
    pythia.readString("TimeShower:QEDshowerByQ = off");
    pythia.readString("ProcessLevel:all = off");
    pythia.readString("HadronLevel:all = on");
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

    std::ostringstream s1, s2, s3, s4;
    s3 << "Random:seed = " << processid;
    s4 << "TimeShower:pTmin = " << Q0;
    std::cout<< s4.str();
    pythia.readString(s3.str());
    pythia.readString(s4.str());
    pythia.init();
}

void PGunWShower::Generate(double pT, std::vector<particle> & plist){
    double x, y;
    TRENToSampler.SampleXY(y, x);
        plist.clear();
        pythia.event.reset();
        pythia.event.append(21, 23, 101, 102, 
                         pT, 0., 0., pT, 0.); 
        pythia.event.append(21, 23, 102, 101, 
                         -pT, 0., 0., pT, 0.); 

        /*pythia.event.append(1, 23, 101, 0, 
                         pT, 0., 0., pT, 0.); 
        pythia.event.append(-1, 23, 0, 101, 
                         -pT, 0., 0., pT, 0.); 
        */
        pythia.event[1].scale(pT);
        pythia.event[2].scale(pT);
        pythia.forceTimeShower(1, 2, pT);
	pythia.next();
        auto & event = pythia.event;
        color_count = event.lastColTag()+1;
	for (size_t i = 0; i < event.size(); ++i) {
		auto p = event[i];
		if (p.isParton() && event[p.daughter1()].isHadron()) {
		//if (p.isFinal()) {
           // final momenta 
            fourvec p0{p.e(), p.px(), p.py(), p.pz()};
            particle _p; 
            _p.pid = p.id();
            _p.mass = std::abs(p.m());
            _p.x0 = fourvec{0,x,y,0};
            _p.x = _p.x0; 
	    double t0 = 0.;
	    reference_pmu(i, t0, event);
            _p.tau_i = t0;
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
            _p.charged = p.isCharged(); 
            plist.push_back(_p);
	    }
	}
}


#endif
