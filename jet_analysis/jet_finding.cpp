#include <fstream>
#include <string>
#include <math.h>
#include <vector>
#include "jet_finding.h"
#include "workflow.h"
#include "lorentz.h"

void FindJet(std::vector<particle> plist, 
             std::vector<current> jlist,
             double jetRadius, 
             double jetpTMin, 
             double jetyMin, 
             double jetyMax,
             std::string fname, double sigma_gen) {
    int power = -1; // power=-1: anti-kT
    //fastjet.init();
    // Define a jet
    fastjet::JetDefinition jetDef(fastjet::genkt_algorithm, jetRadius, power);
    // Objetcts to operate on
    std::vector<fastjet::PseudoJet> fjInputs;
    // cuts
    fastjet::Selector select_rapidity = fastjet::SelectorRapRange(jetyMin, jetyMax);
    for (auto & p : plist){
        fastjet::PseudoJet fp(p.p.x(), p.p.y(), p.p.z(), p.p.t());
        fp.set_user_info(new MyInfo(p.pid, p.origin, sigma_gen));
        fjInputs.push_back(fp);
    }
    
    std::vector<fastjet::PseudoJet> inclusiveJets, jets;
    fastjet::ClusterSequence clustSeq(fjInputs, jetDef);
    inclusiveJets = clustSeq.inclusive_jets(jetpTMin);
    std::vector<fastjet::PseudoJet> selected_jets = select_rapidity(inclusiveJets);
    jets = sorted_by_pt(selected_jets);
    std::ofstream f(fname.c_str(), std::ios_base::app);
    std::vector<fourvec> AllJets;

    // The source list is not part of the jet finding, but used as a 
    // correlated back-ground (may be a confusing name) 
    // to further correct the jet energy.
    double d2 =  std::pow(jetRadius, 2);
    for (auto & j : jets){
        fourvec jetP{j.e(), j.px(), j.py(), j.pz()};
        double pT = jetP.xT();
        double rap = jetP.pseudorap();
        double phi = jetP.phi();

        // The source list is not part of the jet finding, but used as a 
        // correlated back-ground (may be a confusing name) 
        // to further correct the jet energy.
        double cs = 1./std::sqrt(3.);
        int Nphi = 27;
        double dphi = M_PI*2./Nphi;
        double dcostheta = 0.1;
        double dOmega = 1./4./M_PI*dphi*dcostheta;
        for (auto & s : jlist) {
            // correct jet energy if the source 
            // current point inside the jet radius
            double srap = s.x.rap();
            for (int i=0; i<Nphi; i++){
                for (double costheta=-1.; costheta<=1.; costheta+=dcostheta){
                    double iphi = dphi*i;
                    double sintheta = std::sqrt(1.-costheta*costheta);
                    fourvec k{cs, sintheta*std::cos(iphi),
                              sintheta*std::sin(iphi), costheta};
                    double P0 = (k.t()*s.p.t()+k.x()*s.p.x()
                               +k.y()*s.p.y()+k.z()*s.p.z())*dOmega;
                    fourvec dpmu{P0/k.t(), 3*P0*k.x(), 
                                 3*P0*k.y(), 3*P0*k.z()};
                    k = k.boost_back(s.v[0], s.v[1], s.v[2]);
                    double dist2 = std::pow(rap-srap, 2)
                                   +std::pow(phi-k.phi(), 2);
                    if (dist2 < d2){
                        dpmu = dpmu.boost_back(s.v[0], s.v[1], s.v[2]);
                        pT += (dpmu.x()*jetP.x()+dpmu.y()*jetP.y())/(jetP.xT()+1e-9);
                    }
                }
            }
        }  
        f << pT << " " << jetP.phi() << " " 
          << jetP.pseudorap() << " " << jetP.m2() << " " 
          << sigma_gen << " "  
          << jetP << std::endl;
    }   
}

void JetShape(std::vector<particle> plist, 
             std::vector<current> jlist,
             double jetRadius, 
             double jetpTMin, 
             double jetyMin, 
             double jetyMax,
             std::string fname, double sigma_gen) {
    double d2 = jetRadius*jetRadius;
    int power = -1; // power=-1: anti-kT
    //fastjet.init();
    // Define a jet
    fastjet::JetDefinition jetDef(fastjet::genkt_algorithm, jetRadius, power);
    // Objetcts to operate on
    std::vector<fastjet::PseudoJet> fjInputs;
    // cuts
    fastjet::Selector select_rapidity = fastjet::SelectorRapRange(jetyMin, jetyMax);
    for (auto & p : plist){
        fastjet::PseudoJet fp(p.p.x(), p.p.y(), p.p.z(), p.p.t());
        fp.set_user_info(new MyInfo(p.pid, p.origin, sigma_gen));
        fjInputs.push_back(fp);
    }
    
    std::vector<fastjet::PseudoJet> inclusiveJets, jets;
    fastjet::ClusterSequence clustSeq(fjInputs, jetDef);
    inclusiveJets = clustSeq.inclusive_jets(jetpTMin);
    std::vector<fastjet::PseudoJet> selected_jets = select_rapidity(inclusiveJets);
    jets = sorted_by_pt(selected_jets);
    std::ofstream f(fname.c_str(), std::ios_base::app);
    for (auto & j : jets){
        fourvec jetP{j.e(), j.px(), j.py(), j.pz()};
        double pT = jetP.xT();
        double rap = jetP.pseudorap();
        double phi = jetP.phi();

        // The source list is not part of the jet finding, but used as a 
        // correlated back-ground (may be a confusing name) 
        // to further correct the jet energy.
        double cs = 1./std::sqrt(3.);
        int Nphi = 27;
        double dphi = M_PI*2./Nphi;
        double dcostheta = 0.1;
        double dOmega = 1./4./M_PI*dphi*dcostheta;
        for (auto & s : jlist) {
            // correct jet energy if the source 
            // current point inside the jet radius
            double srap = s.x.rap();
            for (int i=0; i<Nphi; i++){
                for (double costheta=-1.; costheta<=1.; costheta+=dcostheta){
                    double iphi = dphi*i;
                    double sintheta = std::sqrt(1.-costheta*costheta);
                    fourvec k{cs, sintheta*std::cos(iphi),
                              sintheta*std::sin(iphi), costheta};
                    double P0 = (k.t()*s.p.t()+k.x()*s.p.x()
                               +k.y()*s.p.y()+k.z()*s.p.z())*dOmega;
                    fourvec dpmu{P0/k.t(), 3*P0*k.x(), 
                                 3*P0*k.y(), 3*P0*k.z()};
                    k = k.boost_back(s.v[0], s.v[1], s.v[2]);
                    double dist2 = std::pow(rap-srap, 2)
                                   +std::pow(phi-k.phi(), 2);
                    if (dist2 < d2){
                        pT += (dpmu.x()*k.x()+dpmu.y()*k.y())/(k.xT()+1e-9);
                    }
                }
            }
        }

        // Jet shapes
        double w1 = 0., w0 = 0.;
        double bins[18] = {0, .05, .1, .15, .2, .25, .3,
                        .35, .4, .45, .5, .6, .7, .8, 1,
                        1.5, 2.0, 3.0};
        if (pT > 120 && std::abs(rap)<1.6){
            std::vector<double> hist, hist0;
            hist.resize(17);
            hist0.resize(17);
            for (auto & it : hist) it = 0.;
            for (auto & it : hist0) it = 0.;

            // particle contribution
            for (auto & p : plist){       
                double jrap = p.p.pseudorap();
                double jphi = p.p.phi();
                double r = std::sqrt(std::pow(phi-jphi,2)
                                    +std::pow(rap-jrap,2));
                if (r < 3){
                    int index = 0;
                    for (index=0; index<17; index++){
                        if (r>bins[index] && r<bins[index+1]) break;
                    }
                    double jpT = p.p.xT();
                    hist[index] += jpT;
                    w1 += jpT;
                    hist0[index] += jpT;
                    w0 += jpT;
                }
            }


            // The source list is not part of the jet finding, but used as a 
            // correlated back-ground (may be a confusing name) 
            // to further correct the jet energy.
            double cs = 1./std::sqrt(3.);
            int Nphi = 27;
            double dphi = M_PI*2./Nphi;
            double dcostheta = 0.1;
            double dOmega = 1./4./M_PI*dphi*dcostheta;
            for (auto & s : jlist) {
                // correct jet energy if the source 
                // current point inside the jet radius
                double srap = s.x.rap();
                for (int i=0; i<Nphi; i++){
                    for (double costheta=-1.; costheta<=1.; costheta+=dcostheta){
                        double iphi = dphi*i;
                        double sintheta = std::sqrt(1.-costheta*costheta);
                        fourvec k{cs, sintheta*std::cos(iphi),
                                  sintheta*std::sin(iphi), costheta};
                        double P0 = (k.t()*s.p.t()+k.x()*s.p.x()
                                  +k.y()*s.p.y()+k.z()*s.p.z())*dOmega;
                        fourvec dpmu{P0/k.t(), 3*P0*k.x(), 
                                  3*P0*k.y(), 3*P0*k.z()};
                        k = k.boost_back(s.v[0], s.v[1], s.v[2]);
                        double dist = std::sqrt(std::pow(rap-srap, 2)
                                               +std::pow(phi-k.phi(), 2));
                        if (dist < 3.){
                            int index=0;
                            for (index=0; index<17; index++){
                                if (dist>bins[index]&&dist<bins[index+1]) 
                                    break;
                            }
                            dpmu = dpmu.boost_back(s.v[0], s.v[1], s.v[2]);
                            double jpT = (dpmu.x()*k.x()+dpmu.y()*k.y())/(k.xT()+1e-9);
                            hist[index] += jpT;
                            w1 += jpT;
                        }
                    }   
                }
            }

            for (auto & it : hist) f << it/w1 << " ";
            for (auto & it : hist0) f << it/w0 << " ";
            f << sigma_gen << std::endl;
        }
    }
}

