#ifndef PYTHIA_WRAPPER_H
#define PYTHIA_WRAPPER_H

#include "simpleLogger.h"
#include "Pythia8/Pythia.h"
#include "workflow.h"
#include <sstream>
#include <unistd.h>
#include "predefine.h"
using namespace Pythia8;

class HQGenerator{
public: 
    HQGenerator(std::string f_pythia, std::string f_trento, 
                int iev, double pTHL, double pTHH, double _Q0);
    void Generate(std::vector<particle> & plist, int Neve, double ycut);
private:
    double pL, pH, Q0, sigma0;
    Pythia pythia;
    std::string f_pythia;
    TransverPositionSampler TRENToSampler;
};

HQGenerator::HQGenerator(std::string f_p, std::string f_trento, int iev, double pTHL, double pTHH, double _Q0):
pL(pTHL), pH(pTHH),f_pythia(f_p),TRENToSampler(f_trento, iev)
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

void reference_pmu(int i, double & tau, Event & event){
    auto p = event[i];
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


void HQGenerator::Generate(std::vector<particle> & plist, int Neve, double ycut){
    double x0, y0;
    int Ncharm = 0, Nbottom = 0;
    plist.clear();
    for(int Ntot=0; Ntot<Neve; Ntot++){
        if (Ntot%1000==0) LOG_INFO << "Ncharm = " << Ncharm << ", Nbottom = " << Nbottom;
        TRENToSampler.SampleXY(y0, x0);
        pythia.next();
        double weight = pythia.info.sigmaGen()/Neve;
        for (size_t i = 0; i < pythia.event.size(); ++i) {
            auto p = pythia.event[i];
            bool triggered = ((p.idAbs() == 5) || (p.idAbs() == 4))
                      && p.isFinal() && (std::abs(p.y())< ycut);
            if (triggered) {
	        double t0 = 0.;
	        reference_pmu(i, t0, pythia.event);
                if (p.idAbs() == 4) Ncharm ++;
                if (p.idAbs() == 5) Nbottom ++;
                for(int iphi=0; iphi<8; iphi++){
                    double phi = iphi*2*M_PI/8;
                    double cos = std::cos(phi), sin = std::sin(phi);
                    fourvec p0{p.e(), p.px()*cos-p.py()*sin, 
                                      p.px()*sin+p.py()*cos, p.pz()};
                    particle _p; 
                    _p.pid = p.idAbs();
                    _p.mass = std::abs(p.m());
                    _p.x0 = fourvec{0,x0,y0,0};
                    _p.tau_i = t0;
                    _p.x = _p.x0; 
                    _p.Q0 = Q0;
                    _p.Q00 = Q0;
                    _p.p0 = p0;
                    _p.col = p.col();
                    _p.acol = p.acol();
                _p.p = _p.p0; 
                _p.weight = weight/8;
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
    }
    //double normW = 0.5*pythia.info.sigmaGen()/pythia.info.weightSum();
    //for (auto & p : plist) p.weight *= normW;
    //LOG_INFO << "norm = " << normW;
}


#endif
