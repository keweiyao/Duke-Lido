// main21.cc is a part of the PYTHIA pythia.event generator.
// Copyright (C) 2018 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program.
// It illustrates how to feed in a single particle (including a resonance)
// or a toy parton-level configurations.

#include "Pythia8/Pythia.h"
#include "Hadronize.h"
#include "lorentz.h"
#include "random.h"
#include "simpleLogger.h"
#include <fstream>
using namespace Pythia8;


JetDenseMediumHadronize::JetDenseMediumHadronize(){
    pythia.readString("ProcessLevel:all = off");
    pythia.readString("TimeShower:QEDshowerByQ = off");
    pythia.readString("Next:numberShowInfo = 0");
    pythia.readString("Next:numberShowProcess = 0");
    pythia.readString("Next:numberShowEvent = 0");
    pythia.readString("1:m0 = 0");
    pythia.readString("2:m0 = 0");
    pythia.readString("3:m0 = 0");
    pythia.readString("4:m0 = 1.3");
    pythia.readString("5:m0 = 4.2");
    pythia.readString("PartonLevel:Remnants = on");
    pythia.readString("HadronLevel:all = on");
    pythia.readString("HadronLevel:Decay = on");
    pythia.readString("StringZ:usePetersonC=on");
    pythia.readString("StringZ:usePetersonB=on");
    pythia.readString("4:mayDecay = off");
    pythia.readString("5:mayDecay = off");
pythia.readString("411:mayDecay = off");
pythia.readString("421:mayDecay = off");
pythia.readString("413:mayDecay = off");
pythia.readString("423:mayDecay = off");
pythia.readString("511:mayDecay = off");
pythia.readString("521:mayDecay = off");
pythia.readString("513:mayDecay = off");
pythia.readString("523:mayDecay = off");
    pythia.init();
}

double GetMass(int pid){
    int absid = std::abs(pid);
    if (absid == 21) return 0.;
    else if (absid == 123 || absid == 1 || absid == 2 || absid == 3) return 0.;
    else if (absid == 4) return 1.3;
    else if (absid == 5) return 4.2;
    else {
        LOG_WARNING << "I don't know what pid = " << pid << " is, but going to give it a zero mass.";
        return 0.;
    }
}

double GetThreshold(int pid){
    int absid = std::abs(pid);
    if (absid == 21) return 0.;
    else if (absid == 123 || absid == 1 || absid == 2 || absid == 3) return 0.;
    else if (absid == 4) return 1.3;
    else if (absid == 5) return 4.2;
    else {
        LOG_WARNING << "I don't know what pid = " << pid << " is, but going to give it a zero mass.";
        return 0.;
    }
}

void FormChain(particle pi, particle pf, 
               std::vector<particle> & chain, 
               std::vector<particle> & plist, 
               std::vector<particle> & thermal){
    int size0 = chain.size();
    for (std::vector<particle>::iterator it = plist.begin(); 
         it != plist.end();){
        if (it->col==pf.acol && pf.acol!=0){
            chain.push_back(*it);
            pf = chain.back();
            it = plist.erase(it);
            continue;
        }
        if (it->acol==pi.col && pi.col!=0){
            chain.push_back(*it);
            pi = chain.back();
            it = plist.erase(it);
            continue;
        }
        else it++;
    }
    int size1 = chain.size();
    if (size1==size0) {
        // fill open end point with thermal particles
        if (pi.col!=0){
            particle th;
            th.col = 0;
            th.acol = pi.col;
	    double vx = 0.;
	    double vy = 0.;
	    double vz = pi.x.z()/(1e-9+pi.x.t());
            if (std::abs(vz)>.9999) vz = vz/std::abs(vz)*.9999;
            th.pid = -std::abs(Srandom::sample_flavor(3));
	    th.p = Srandom::generate_thermal_parton_with_boost(
                  std::max(pi.Tf,.165), vx, vy, vz);
            th.mass = 0.;
	    th.x0 = pi.x;
            th.vcell.resize(3);
	    th.vcell[0] = vx;
            th.vcell[1] = vy;
            th.vcell[2] = vz;
	    th.x = pi.x;
	    th.Tf = pi.Tf;
            th.Q0 = 0.;
            thermal.push_back(th);
            chain.push_back(th);
        }
        // fill open end point with thermal particles
        if (pf.acol!=0){
            particle th;
            th.col = pf.acol;
            th.acol = 0;
	    double vx = 0.;
            double vy = 0.;
            double vz = pf.x.z()/(1e-9+pf.x.t());
	    if (std::abs(vz)>.9999) vz = vz/std::abs(vz)*.9999;
            th.pid = std::abs(Srandom::sample_flavor(3));
	    th.p = Srandom::generate_thermal_parton_with_boost(
                  std::max(pf.Tf,.165), vx, vy, vz);
            th.mass = 0.;
	    th.x0 = pf.x;
            th.vcell.resize(3);
	    th.vcell[0] = vx;
            th.vcell[1] = vy;
            th.vcell[2] = vz;
	    th.x = pf.x;
	    th.Tf = pf.Tf;
            th.Q0 = 0.;
            thermal.push_back(th);
            chain.push_back(th);
        }
        return;
    }
    else {
        FormChain(pi, pf, chain, plist, thermal);
    }
}

int JetDenseMediumHadronize::hadronize(std::vector<particle> partons, 
                                       std::vector<particle> & hadrons, 
                                       std::vector<particle> & thermal_partons,
                                       double Q0,
                                       int level){
    int Npartons = partons.size();
    const int status = 23;
    hadrons.clear();
    thermal_partons.clear();
    int Ninit, Nfinal;
    // step1: group partons into patches connected by color index
    // step2: for each color chain, 
    //        sample thermal partons to form color singlet
    // step3: call pythia force-time-shower with hadronization
    // define a list of color chains
    std::vector<std::vector<particle> > chains;
    while(!partons.empty()){
        std::vector<particle> chain;
        chain.clear();
        chain.push_back(partons.back());
        partons.pop_back();
        FormChain(chain.front(), chain.back(), chain, 
                  partons, thermal_partons);
        chains.push_back(chain);
    }
   
    for (auto & c : chains){
        double maxQ0 = Q0;
        pythia.event.reset();
        int count=1;
        for (auto &p : c){
             pythia.event.append(p.pid, 23, p.col, p.acol, 
		                 p.p.x(), p.p.y(), p.p.z(), p.p.t(), p.mass);
             pythia.event[count].scale(p.Q0);
             count++;
             maxQ0 = (p.Q0 > maxQ0) ? p.Q0 : maxQ0;
        }     
        pythia.forceTimeShower(1,count-1,maxQ0);
	pythia.next();
        int Nff=0;
        for (int i = 0; i < pythia.event.size(); ++i) {
            auto ip = pythia.event[i];
            bool good = false;
            if (level==1){
                good = ip.isFinal();
            }
            if (level==0){
                good = (ip.isParton() 
                        && pythia.event[ip.daughter1()].isHadron())
                      || (ip.isFinal() && ip.isParton());
            }
            if (good) {
                Nff ++;
                particle h;
                h.pid = ip.id();
                h.p.a[0] = ip.e();
                h.p.a[1] = ip.px();
                h.p.a[2] = ip.py();
                h.p.a[3] = ip.pz();
                h.mass = std::sqrt(h.p.t()*h.p.t() - h.p.pabs2());
                h.weight = 1;
		h.charged = ip.isCharged();
                hadrons.push_back(h);
            }
            
            if (level==1 && ip.isFinal() && ip.isParton()){
                //LOG_INFO << "recombin needed" << ip.id() << " " 
		//	 << ip.e() << " " << ip.px() << " "
		//	 << ip.py() << " " << ip.pz();
            }
        }
    }        
    //LOG_INFO << Npartons << " to " << hadrons.size();
    return hadrons.size();
}
