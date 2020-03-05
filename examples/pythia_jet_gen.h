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
              int pTHL, int pTHH, int iev);
    void Generate(std::vector<particle> & plist, int heavyid);
    double sigma_gen(void){
        return pythia.info.sigmaGen();
    }
private:
    Pythia pythia;
    TransverPositionSampler TRENToSampler;
};

PythiaGen::PythiaGen(std::string f_pythia, std::string f_trento,
                     int pTHL, int pTHH, int iev):
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

void find_production_x(int i, fourvec & x, Event & event){
    double tf;
    auto p = event[i];
    int im1 = p.mother1();
    int im2 = p.mother2();
    if (im1==im2 && im2 > 0){
        // simply recoil effect, go on
        find_production_x(im1, x, event);
    }
    if (im1==0 && im2 ==0) {
        return;
    }
    if (im1>0 && im2==0){
        auto m1 = event[im1];
        // FSR or ISR

        if (p.status() == 51) {
             int d1 = m1.daughter1();
             int d2 = m1.daughter2();
             auto D1 = event[d1];
             auto D2 = event[d2];
             
             double ME = D1.e()+D2.e();
             double f = p.e()/ME;
             double Mpx = D1.px()+D2.px(), Mpy = D1.py()+D2.py(), Mpz = D1.pz()+D2.pz();
             double p2 = p.px()*p.px() + p.py()*p.py() + p.pz()*p.pz();
             double Q2 = Mpx*Mpx + Mpy*Mpy + Mpz*Mpz;
             double pdotQ = p.px()*Mpx + p.py()*Mpy + p.pz()*Mpz;
             double pT2 = p2 - pdotQ*pdotQ/Q2;
             tf = 2.*f*(1.-f)*ME/pT2;
       
            x.a[0] = x.t() + tf;
            x.a[1] = x.x() + tf*Mpx/ME;
            x.a[2] = x.y() + tf*Mpy/ME;
            x.a[3] = x.z() + tf*Mpz/ME;
            find_production_x(im1, x, event);
        }
        else{
            find_production_x(im1, x, event);
        }
    }
    if (im1 != im2 && im1 > 0 && im2 > 0){
        // hard 2->n process, time scale is short, and we negelect
        return;
    }
}

void PythiaGen::Generate(std::vector<particle> & plist, int heavyid){
    double x, y;
    TRENToSampler.SampleXY(x, y);
    LOG_INFO << x/5.076 << " " << y/5.076; 
    bool triggered = false;
    do{
        if (heavyid < 4) triggered = true;
        plist.clear();
        pythia.next();
        double weight = pythia.info.sigmaGen();
	for (size_t i = 0; i < pythia.event.size(); ++i) {
		auto p = pythia.event[i];
                
		if (p.isFinal() && std::abs(p.y())< 3.5) {
		    // final momenta 
		    fourvec p0{p.e(), p.px(), p.py(), p.pz()};
		    particle _p; 
		    _p.pid = p.id();
		    _p.mass = std::abs(p.m());
		    _p.x0 = fourvec{0,x,y,0};
                    _p.x = _p.x0; 
		    //find_production_x(i, _p.x, pythia.event); 
		    _p.p0 = p0;
 
		    if (std::abs(_p.pid) != 4 && std::abs(_p.pid) != 5) {
		        _p.mass = 0;
		    }
                    if (std::abs(_p.pid) == heavyid) triggered = true;
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
