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
#include <algorithm>
using namespace Pythia8;


JetDenseMediumHadronize::JetDenseMediumHadronize(){
    pythia.readString("Tune:pp=19");
    pythia.readString("ProcessLevel:all = off");
    pythia.readString("Print:quiet = on");
    pythia.readString("SoftQCD:all = off");
    pythia.readString("PromptPhoton:all=off");
    pythia.readString("WeakSingleBoson:all=off");
    pythia.readString("WeakDoubleBoson:all=off");
    pythia.readString("SpaceShower:QEDshowerByQ=off");
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
    pythia.readString("HadronLevel:Decay = off");
    pythia.readString("StringZ:usePetersonC=on");
    pythia.readString("StringZ:usePetersonB=on");
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
	       particle & front, particle & back,
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
	front = pi;
	back = pf;
        return;
    }
    else FormChain(pi, pf, front, back, chain, plist, thermal);
}

int JetDenseMediumHadronize::hadronize(std::vector<particle> partons, 
                                       std::vector<particle> & hadrons, 
                                       std::vector<particle> & thermal_partons,
                                       double Q0, double Tf,
                                       int level, int Nos){
    std::vector<particle> colorless_ensemble;
    for(auto & p : partons) colorless_ensemble.push_back(p);
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
    std::vector<particle> endpoint_cols, endpoint_acols;
    endpoint_cols.clear();
    endpoint_acols.clear();
    while(!partons.empty()){
	particle front, back;
        std::vector<particle> chain;
        chain.clear();
        chain.push_back(partons.back());
        partons.pop_back();
        FormChain(chain.front(), chain.back(), 
		  front, back,
		  chain, 
                  partons, thermal_partons);
        chains.push_back(chain);
	endpoint_cols.push_back(front);
	endpoint_acols.push_back(back);
    }
    // Sort the endpoints by rapidty
    std::sort(endpoint_cols.begin(), endpoint_cols.end(), 
      [](const particle & a, const particle & b){
          return a.x.rap()>b.x.rap();
      }
      );
    std::sort(endpoint_acols.begin(), endpoint_acols.end(), 
      [](const particle & a, const particle & b){
          return a.x.rap()>b.x.rap();
      }
      );

    // find out the endpoint that are not color neutral
    // connect them randomly to thermal partons
    for (int i=0; i<endpoint_cols.size();i++){
	auto & p1 = endpoint_cols[i];
	auto & p2 = endpoint_acols[i];
	
	int acol = p1.col;
	int col = p2.acol;
	int pid;
	double tau, rap;
	    if (col==0 && acol==0) continue;
	    if (col==0 && acol!=0){
                pid = -std::abs(Srandom::sample_flavor(3));
		tau = p1.x.tau();
                rap = p1.x.rap();
	    }
	    if (col!=0 && acol==0){
                pid = std::abs(Srandom::sample_flavor(3));
                tau = p2.x.tau();
                rap = p2.x.rap();
	    }
            if (col!=0 && acol!=0){
                pid = 21;
		tau = (p1.x.tau()+p2.x.tau())/2.;
                rap = (p1.x.rap()+p2.x.rap())/2.;
	    }
            particle th;
	    th.pid = pid;
            th.col = col;
            th.acol = acol;
            th.x0 = fourvec{tau*std::cosh(rap),0.,0.,tau*std::sinh(rap)};
	    th.x = th.x0;
	    double vzgrid = th.x.z()/th.x.t();
	    do{
            th.p = Srandom::generate_thermal_parton_with_boost(
                           Tf, 0, 0, vzgrid);
	    }while(th.p.xT()>.5);
            th.mass = 0.;
            th.vcell.resize(3);
            th.vcell[0] = 0.;
            th.vcell[1] = 0.;
            th.vcell[2] = vzgrid;
            th.Tf = Tf;
            th.Q0 = 0.5;
            thermal_partons.push_back(th);
	    colorless_ensemble.push_back(th);
    }
    for (int i=0; i<Nos; i++){
    double maxQ0 = Q0;
    pythia.event.reset();
    int count=1;
    for (auto & p : colorless_ensemble){
         pythia.event.append(p.pid, 23, p.col, p.acol, 
			 p.p.x(), p.p.y(), p.p.z(), p.p.t(), p.mass);
         pythia.event[count].scale(p.Q0);
	 maxQ0 = (p.Q0>maxQ0) ? p.Q0 : maxQ0;
         count++;
    }     
    pythia.forceTimeShower(1,count-1,maxQ0);
    pythia.next();
    int Nff=0;
    for (int i = 0; i < pythia.event.size(); ++i) {
        auto ip = pythia.event[i];
        bool good = false;
	int absid = std::abs(ip.id());
        if (level==1) good = ip.isFinal();
        if (level==0) good = (ip.isParton() 
                        && pythia.event[ip.daughter1()].isHadron())
                      || (ip.isFinal() && ip.isParton());
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
            
       } 
    }
    return hadrons.size();
}
