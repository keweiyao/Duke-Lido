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
    pythia.readString("HadronLevel:Decay = off");
    /*pythia.readString("4:mayDecay = off");
    pythia.readString("5:mayDecay = off");
    pythia.readString("111:mayDecay = off");
    pythia.readString("211:mayDecay = off");
    pythia.readString("311:mayDecay = off");
    pythia.readString("321:mayDecay = off");
    pythia.readString("411:mayDecay = off");
    pythia.readString("421:mayDecay = off");
    pythia.readString("413:mayDecay = off");
    pythia.readString("423:mayDecay = off");
    pythia.readString("511:mayDecay = off");
    pythia.readString("521:mayDecay = off");
    pythia.readString("513:mayDecay = off");
    pythia.readString("523:mayDecay = off");*/
    pythia.init();
}

fourvec generate_thermal_parton_with_boost(double T, double vx, double vy, double vz){
    // randomly sample the four momentum
    double se = T*Srandom::sample_E_over_T(Srandom::gen);
    double phi = Srandom::dist_phi(Srandom::gen);
    double cos = Srandom::dist_costheta(Srandom::gen);
    double sin = std::sqrt(1.-cos*cos);
    fourvec p{se, se*sin*std::cos(phi), se*sin*std::sin(phi),se*cos};
    return p.boost_back(vx, vy, vz);
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

int JetDenseMediumHadronize::hadronize(std::vector<particle> partons, 
                                       std::vector<particle> & hadrons, 
                                       std::vector<particle> & thermal_partons,
                                       double Q0,
                                       int level){
    
    const int color1 = 101, color2 = 102;
    const int status = 23;
    hadrons.clear();
    thermal_partons.clear();
    for (auto & p : partons){
        std::vector<particle> plist;
        plist.push_back(p);
        pythia.event.reset();
        // velocity of the comoving frame of the medium cell 
        // where parton pass throught the QGP phase boundary (not very well defined hyper-surface)
        double vx = p.vcell[0];
        double vy = p.vcell[1];
        double vz = p.vcell[2];
        
        double weight = p.weight;
        double T = std::max(p.Tf, 0.1);
        //LOG_INFO << T << " " << vx << " " << vy << " " << vz;
        int absid = std::abs(p.pid);
 
        
        // quark case
        double threshold = GetThreshold(p.pid), sqrts;
        do{
        fourvec ptot{0., 0., 0., 0.};
        if (absid==123 || absid==1 || absid==2 
           || absid==3 || absid==4 || absid==5) {
            int cpid, color, anticolor;
            if (p.pid>0) {
                // a quark
                color = color1;
                anticolor = 0;
                cpid = -(1+std::rand()%3); // randomly choose its flavor from u d s
            }
            else{
                // an anti quark
                color = 0;
                anticolor = color2;
                cpid = (1+std::rand()%3); // randomly choose its flavor from u d s
            }
            // attach the hard parton
            auto pth = generate_thermal_parton_with_boost(T, vx, vy, vz);
            pythia.event.append(p.pid, status, color, anticolor, 
                         p.p.x(), p.p.y(), p.p.z(), p.p.t(), p.mass); 
            pythia.event.append(cpid, status, anticolor, color, 
                         pth.x(), pth.y(), pth.z(), pth.t(), 0.0); 
	    pythia.event[1].scale(p.Q0);
            pythia.event[2].scale(0.4);
            pythia.forceTimeShower(1,2,p.Q0);
            particle thermal_p;
            thermal_p.pid = cpid;
            thermal_p.p = pth;
            thermal_p.mass = 0.;
            thermal_p.x0 = p.x;
            thermal_p.vcell = p.vcell;
            thermal_p.x = p.x;
            thermal_p.Tf = T;
            thermal_partons.push_back(thermal_p);
            ptot = thermal_p.p + p.p;
            plist.push_back(thermal_p);
        }
        if (absid == 21) {
            int color = color1;
            int anticolor = color2;
            int pid1 = (1+std::rand()%3);
            int pid2 = -(1+std::rand()%3);
            auto pth1 = generate_thermal_parton_with_boost(T, vx, vy, vz);
            auto pth2 = generate_thermal_parton_with_boost(T, vx, vy, vz);
            pythia.event.append(p.pid, status, color, anticolor, 
                         p.p.x(), p.p.y(), p.p.z(), p.p.t(), 0.0); 
            pythia.event.append(pid1, status, anticolor, 0, 
                         pth1.x(), pth1.y(), pth1.z(), pth1.t(), 0.0);  
            pythia.event.append(pid2, status, 0, color, 
                         pth2.x(), pth2.y(), pth2.z(), pth2.t(), 0.0);  
            pythia.event[1].scale(p.Q0);
	    pythia.event[2].scale(0.4);
	    pythia.event[3].scale(0.4);
	    pythia.forceTimeShower(1,3,p.Q0);
            particle thermal_p1;
            thermal_p1.pid = pid1;
            thermal_p1.p = pth1;
            thermal_p1.mass = 0.;
            thermal_p1.x0 = p.x;
            thermal_p1.vcell = p.vcell;
            thermal_p1.x = p.x;
            thermal_p1.Tf = T;
            thermal_partons.push_back(thermal_p1);
            plist.push_back(thermal_p1);

            particle thermal_p2;
            thermal_p2.pid = pid2;
            thermal_p2.p = pth2;
            thermal_p2.mass = 0.;
            thermal_p2.x0 = p.x;
            thermal_p2.vcell = p.vcell;
            thermal_p2.x = p.x;
            thermal_p2.Tf = T;
            thermal_partons.push_back(thermal_p2);
            plist.push_back(thermal_p2);

            ptot = thermal_p1.p + thermal_p2.p + p.p;
        }
        sqrts = std::sqrt(dot(ptot, ptot));
        //LOG_INFO << "sqrts of fragmentation = " << sqrts << " GeV"; 
        }while(sqrts < threshold);

        int s = pythia.next();
        bool need_recombiantion = false;
        for (int i = 0; i < pythia.event.size(); ++i) {
            auto ip = pythia.event[i];
            bool good = false;
            if (level==1){
                good = ip.isFinal() && ip.isHadron();
            }
            if (level==0){
                good = (ip.isParton() 
                       && pythia.event[ip.daughter1()].isHadron())
                      || (ip.isFinal() && ip.isParton());
            }
            if (good) {
                particle h;
                h.pid = ip.id();
                h.p.a[0] = ip.e();
                h.p.a[1] = ip.px();
                h.p.a[2] = ip.py();
                h.p.a[3] = ip.pz();
                h.mass = std::sqrt(h.p.t()*h.p.t() - h.p.pabs2());
                h.weight = 1;
                h.x0 = p.x;
                h.x = p.x;
                h.vcell = p.vcell;
                hadrons.push_back(h);
            }
            
            if (level==1 && ip.isFinal() && ip.isParton()){
                need_recombiantion = true;
                break;
            }
        }
        /*if (level==1 && need_recombiantion){
            // go back to plist:
            fourvec ptot{0.,0.,0.,0.};
            for (auto & ip : plist) ptot = ptot + ip.p;
            particle h;
            h.mass = -1;
            h.p.a[1] = ptot.x();
            h.p.a[2] = ptot.y();
            h.p.a[3] = ptot.z();

            if (plist[0].pid==21) plist.erase(plist.begin());

            if (std::abs(plist[0].pid) == 5) {
                if (std::abs(plist[1].pid) == 1){
                    h.mass = 5.27929;
                    h.pid = (plist[0].pid>0)?(-521):(521);
                }
                if (std::abs(plist[1].pid) == 2){
                    h.mass = 5.27961;
                    h.pid = (plist[0].pid>0)?(-511):(511);
                }
                if (std::abs(plist[1].pid) == 3){
                    h.mass = 5.36679;
                    h.pid = (plist[0].pid>0)?(-531):(531);
                }
            }
            else if (std::abs(plist[0].pid) == 4) {
                if (std::abs(plist[1].pid) == 1){
                    h.mass = 1.86962;
                    h.pid = (plist[0].pid>0)?(421):(-421);
                }
                if (std::abs(plist[1].pid) == 2){
                    h.mass = 1.864874;
                    h.pid = (plist[0].pid>0)?(411):(-411);
                }
                if (std::abs(plist[1].pid) == 3){
                    h.mass = 1.96847;
                    h.pid = (plist[0].pid>0)?(431):(-431);
                }
            }
            else if (std::abs(plist[0].pid) == 3) {
                if (std::abs(plist[1].pid) == 1){
                    h.mass = 0.493677;
                    h.pid = (plist[0].pid>0)?(-321):(321);
                }
                if (std::abs(plist[1].pid) == 2){
                    h.mass = 0.497611;
                    h.pid = (plist[0].pid>0)?(-311):(311);
                }
                if (std::abs(plist[1].pid) == 3){
                    h.mass = 1.019;
                    h.pid = (plist[0].pid>0)?(-331):(331);
                }
            }
            h.p.a[0] = std::sqrt(h.mass*h.mass + h.p.pabs2());
            h.weight = -1;
            h.x0 = p.x;
            h.x = p.x;
            h.vcell = p.vcell;
            hadrons.push_back(h);
        }*/
        
    }
    return hadrons.size();
}
