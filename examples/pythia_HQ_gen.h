#ifndef PYTHIA_WRAPPER_H
#define PYTHIA_WRAPPER_H

#include "simpleLogger.h"
#include "Pythia8/Pythia.h"
#include "workflow.h"
#include <sstream>
#include <unistd.h>
using namespace Pythia8;

class HQGenerator{
public: 
    HQGenerator(std::string f_pythia, std::string f_trento, 
                int iev, double pTHL, double pTHH);
    void Generate(std::vector<particle> & plist, int Neve, double ycut);
private:
    double pL, pH;
    Pythia pythia;
    std::string f_pythia;
    TransverPositionSampler TRENToSampler;
};

HQGenerator::HQGenerator(std::string f_p, std::string f_trento, int iev, double pTHL, double pTHH):
pL(pTHL), pH(pTHH),f_pythia(f_p),TRENToSampler(f_trento, iev)
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
    int processid = getpid();

    std::ostringstream s1, s2, s3;
    s1 << "PhaseSpace:pTHatMin = " << pTHL;
    s2 << "PhaseSpace:pTHatMax = " << pTHH;
    s3 << "Random:seed = " << processid;
    pythia.readString(s1.str());
    pythia.readString(s2.str());
    pythia.readString(s3.str());
    // Init
    pythia.init(); 
    for (int i=0; i<1000; i++) pythia.next();
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

void HQGenerator::Generate(std::vector<particle> & plist, int Neve, double ycut){
    double x0, y0;
    int Ncharm = 0, Nbottom = 0;
    plist.clear();
    for(int Ntot=0; Ntot<Neve; Ntot++){
        if (Ntot%1000==0) LOG_INFO << "Ncharm = " << Ncharm << ", Nbottom = " << Nbottom;
        TRENToSampler.SampleXY(x0, y0);
        pythia.next();
        double weight = pythia.info.sigmaGen()/Neve;
        for (size_t i = 0; i < pythia.event.size(); ++i) {
            auto p = pythia.event[i];
            bool triggered = ((p.idAbs() == 5) || (p.idAbs() == 4))
                      && p.isFinal() && (std::abs(p.y())< ycut);
            if (triggered) {
                if (p.idAbs() == 4) Ncharm ++;
                if (p.idAbs() == 5) Nbottom ++;
                for(int iphi=0; iphi<8; iphi++){
                    double phi = iphi*2*M_PI/8.;
                    double cos = std::cos(phi), sin = std::sin(phi);
                    fourvec p0{p.e(), p.px()*cos-p.py()*sin, 
                                      p.px()*sin+p.py()*cos, p.pz()};
                    particle _p; 
                    _p.pid = p.idAbs();
                    _p.mass = std::abs(p.m());
                    _p.x0 = fourvec{0,x0,y0,0};
                    fourvec xx{0,0,0,0};
                    find_production_x(i, xx, pythia.event);
                    _p.tau_i = xx.t();
                    _p.x = _p.x0; 
                    _p.p0 = p0;
                _p.p = _p.p0; 
                _p.weight = weight/8.;
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
        }
    }
    //double normW = 0.5*pythia.info.sigmaGen()/pythia.info.weightSum();
    //for (auto & p : plist) p.weight *= normW;
    //LOG_INFO << "norm = " << normW;
}


#endif
