#ifndef PYTHIA_WRAPPER_H
#define PYTHIA_WRAPPER_H

#include "simpleLogger.h"
#include "Pythia8/Pythia.h"
#include "workflow.h"
#include <sstream>
#include "predefine.h"
#include "random.h"

using namespace Pythia8;

class HQGenerator{
public: 
    HQGenerator(std::string f_pythia, std::string f_trento, 
              double pTHL, double pTHH, int iev, double _Q0);
    void Generate(std::vector<particle> & plist, int Neve, double ycut);
    double sigma_gen(void){
        return sigma0;
    }
    double maxPT(void){
	return pythia.event.scale();
    }
    fourvec x0(void){
	return _x0;
    }
private:
    double pL, pH, Q0, sigma0;
    Pythia pythia;
    TransverPositionSampler TRENToSampler;
    fourvec _x0;
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

HQGenerator::HQGenerator(std::string f_pythia, std::string f_trento, 
              double pTHL, double pTHH, int iev, double _Q0):
TRENToSampler(f_trento, iev)
{   
    Q0 = _Q0;
    // read pythia settings
    pythia.readFile(f_pythia);
    // suppress output
    pythia.readString("Print:quiet = on");
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
    for (int i=0; i<100; i++) pythia.next();
        sigma0 = pythia.info.sigmaGen();
}

void HQGenerator::Generate(std::vector<particle> & plist, int Neve, double ycut){
    double x, y;
    int Ncharm = 0, Nbottom = 0;
    TRENToSampler.SampleXY(y, x);
    plist.clear();
    for(int Ntot=0; Ntot<Neve; Ntot++){
        if (Ntot%1000==0) LOG_INFO << "Ncharm = " << Ncharm 
                                   << ", Nbottom = " << Nbottom;
        TRENToSampler.SampleXY(y, x);
        pythia.next();
        double weight = pythia.info.sigmaGen()/Neve;
        auto & event = pythia.event;
        for (size_t i = 0; i < event.size(); ++i) {
            auto p = event[i];
            int absid = p.idAbs();
            bool triggered = ((absid == 5) || (absid == 4))
                      && p.isFinal() && (std::abs(p.y())< ycut);
            if (triggered) {
                if (absid == 4) Ncharm ++;
                if (absid == 5) Nbottom ++;
                // Get formation time in Lab frame
                double tForm = 0.;
                reference_pmu(i, tForm, event);
                for (int iphi=0; iphi<12; iphi++){
                    particle _p; 
                    double phi = iphi*2.*M_PI/12.;
                    double cos = std::cos(phi), sin = std::sin(phi);
                    // final momenta 
                    fourvec p0{p.e(), p.px()*cos-p.py()*sin, 
                                      p.px()*sin+p.py()*cos, p.pz()};
                    // assign pid
                    _p.pid = absid;
                    // copy charge
                    _p.charged = p.isCharged(); 
                    // assign mass
                    _p.mass = pid2mass(_p.pid);
                    // put on shell
                    p0 = put_on_shell(p0, _p.pid);
                    // space-time rap:
                    double etas = .5*std::log((p0.t()+p0.z())/(p0.t()-p0.z()));
                    // space-time coordinate
                    _p.x0 = coordinate{0., x, y, etas};
                    _p.x = _p.x0;
                    // initialize momentum in coordinate co-moving frame
                    _p.p0 = p0.boost_to(0., 0., p0.z()/p0.t());
                    _p.p = _p.p0; 
                    // Transform to foramtion time in proper time
                    _p.tau0 = tForm*_p.p.t()/p0.t();
                    // Virtuality of particle production
                    _p.Q0 = Q0;
                    _p.Q00 = Q0;
                    // Copy color indices
                    _p.col = p.col();
                    _p.acol = p.acol();
                    // These are real particles to the transport eqautions
                    _p.is_virtual = false;
                    // Temperature of medium at production
                    _p.T0 = 0.;
                    // Mean-free path
                    _p.mfp0 = 0.;
                    // Velocity of medium
                    _p.vcell.resize(3);
                    _p.vcell[0] = 0.; 
                    _p.vcell[1] = 0.; 
                    _p.vcell[2] = 0.; 
                    // clear its radiation list
                    _p.radlist.clear();
                    // ready to go
                    _p.weight = weight/12.;
                    plist.push_back(_p);
                }
            }
        }
    }
}


#endif
