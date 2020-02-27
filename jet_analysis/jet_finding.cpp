#include <fstream>
#include <string>
#include <math.h>
#include <vector>
#include "jet_finding.h"
#include "workflow.h"
#include "lorentz.h"
#include "simpleLogger.h"

inline int corp_index(double x, double xL, double xH, double dx, int Nx){
    if (x<xL || x>xH) return -1;
    else if (x==xH) return Nx-1;
    else return int((x-xL)/dx);
}

void FindJetTower(std::vector<particle> plist, 
             std::vector<current> jlist,
             std::vector<double> Rs,
             double jetpTMin, 
             double jetyMin, 
             double jetyMax,
             std::string fname, double sigma_gen) {
    // contruct four momentum tower in the eta-phi plane
    int Neta = 200, Nphi = 200;
    double etamin = -3., etamax = 3.;
    double deta = (etamax-etamin)/Neta; 
    // Phi range has to be the same as the range of std::atan2(*,*);
    double phimin = -M_PI, phimax = M_PI; 
    double dphi = (phimax-phimin)/Nphi;
    double dOmega = 1./4./M_PI*dphi;
    std::vector<std::vector<fourvec> > towers;
    towers.resize(Neta);
    for (auto & it : towers) {
        it.resize(Nphi);
        for (auto & iit: it) iit = fourvec{0.,0.,0.,0.};
    }
    // put hard particles into the tower
    for (auto & p : plist){
        double eta = p.p.pseudorap();
        double phi = p.p.phi();
        int ieta = corp_index(eta, etamin, etamax, deta, Neta);
        if (ieta<0) continue;
        int iphi = corp_index(phi, phimin, phimax, dphi, Nphi);
        towers[ieta][iphi] = towers[ieta][iphi] + p.p;
    }
    // source integral contribution to each tower

    for (auto & s : jlist) {
        int ieta = corp_index(s.x.rap(), etamin, etamax, deta, Neta);
        if (ieta<0) continue;
        for (int iphi=0; iphi<Nphi; iphi++){
            double phi = phimin+iphi*dphi;
            double cphi = std::cos(phi), sphi = std::sin(phi);
            double P0 = ( M_PI/2.*s.cs*s.p.t()
                        + cphi*4./3.*s.p.x()
                        + sphi*4./3.*s.p.y()
                        )*dOmega; 
            fourvec k{1., s.cs*cphi, s.cs*sphi, 0};   
            double newphi = k.boost_back(s.v[0], s.v[1], 0).phi();   
            int newiphi = corp_index(newphi, phimin, phimax, dphi, Nphi);
            fourvec dpmu{P0/s.cs, 3*P0*cphi, 3*P0*sphi, 0};
            towers[ieta][newiphi] = towers[ieta][newiphi]+
                                dpmu.boost_back(s.v[0], s.v[1], 0).boost_back(0,0, s.v[2]);
        }
    }
    // Use the towers to do jet finding: anti-kT
    int power = -1; 
    // Define a jet
    for (int iR=0; iR<Rs.size(); iR++){

        double jetRadius = Rs[iR];
        fastjet::JetDefinition jetDef(
                    fastjet::genkt_algorithm, jetRadius, power);
        std::vector<fastjet::PseudoJet> fjInputs;
        fastjet::Selector select_rapidity = 
              fastjet::SelectorRapRange(jetyMin, jetyMax);
        for (int ieta=0; ieta<Neta; ieta++) {
            double eta = etamin+ieta*deta;
            double kz = std::sinh(eta), k0 = std::cosh(eta);
            for (int iphi=0; iphi<Nphi; iphi++) {
                auto iit = towers[ieta][iphi];
                double phi = phimin+iphi*dphi;
                double kx = std::cos(phi), ky = std::sin(phi);
                double pT = iit.x()*kx + iit.y()*ky;
                // when we do the first jet finding,
                // only use towers with positive contribution
                if (iit.t()>0 && pT>0){
                    fastjet::PseudoJet fp(iit.x(), iit.y(), iit.z(), iit.t());
                    fp.set_user_info(new MyInfo(0, 0, sigma_gen));
                    fjInputs.push_back(fp);
                }
            }
        }
        
        fastjet::ClusterSequence clustSeq(fjInputs, jetDef);
        std::vector<fastjet::PseudoJet> AllJets = 
                                clustSeq.inclusive_jets(jetpTMin);
        std::vector<fastjet::PseudoJet> jets =
                        fastjet::sorted_by_pt(select_rapidity(AllJets));

        // once the jet is find, 
        // recombine all the bins with negative contribution
        std::vector<fourvec> results;
        for (auto & j : jets){
            fourvec jetP{j.e(), j.px(), j.py(), j.pz()};
            fourvec newP{0., 0., 0., 0.};
            double Jphi = jetP.phi();
            double Jeta = jetP.pseudorap();
            for (double eta=Jeta-jetRadius; eta<=Jeta+jetRadius; eta+=deta){
                int ieta = corp_index(eta, etamin, etamax, deta, Neta);
                if (ieta<0) continue;
                double Rp = std::sqrt(1e-9+std::pow(jetRadius,2)
                                          -std::pow(Jeta-eta,2)); 

                for (double phi=Jphi-Rp; phi<=Jphi+Rp; phi+=dphi){
                    double newphi = phi;
                    if (phi<-M_PI) newphi+=2*M_PI;
                    if (phi>M_PI) newphi-=2*M_PI;
                    int iphi = corp_index(newphi, phimin, phimax, dphi, Nphi);
                    newP = newP + towers[ieta][iphi];
                }  
            }
            results.push_back(newP);
        }
        std::ofstream f(fname.c_str(), std::ios_base::app);
        for (int j=0; j<results.size(); j++){
            auto jetP = results[j];                  
            f << iR << " " << jetP.xT() << " " << jetP.phi() << " " 
              << jetP.pseudorap() << " " << jetP.m2() << " " 
              << sigma_gen << std::endl;
        }
    }   
}

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

    // The source list is not part of the jet finding, but used as a 
    // correlated back-ground (may be a confusing name) 
    // to further correct the jet energy.
    double d2 =  std::pow(jetRadius, 2);
    std::vector<fourvec> results;
    std::vector<double> corrected_pT;
    for (auto & j : jets){
        fourvec jetP{j.e(), j.px(), j.py(), j.pz()};
        results.push_back(jetP);
        corrected_pT.push_back(jetP.xT());
        

        double pT = jetP.xT();
        double rap = jetP.pseudorap();
        double phi = jetP.phi();
        double cosjet = std::cos(phi);
        double sinjet = std::sin(phi);

    // The source list is not part of the jet finding, but used as a 
    // correlated back-ground (may be a confusing name) 
    // to further correct the jet energy.
    double cs = 1/1.7;
    int Nphi = 10;
    double dphi = 2.*jetRadius/Nphi;
    double dOmega = 1./4./M_PI*dphi;  
    for (auto & s : jlist) {
        // correct jet energy if the source 
        // current point inside the jet radius
        double DeltaRap2 = std::pow(rap - s.x.rap(),2);
        if (DeltaRap2 >= d2) continue;
        for (double iphi=phi-jetRadius; 
            iphi<=phi+jetRadius; iphi+=dphi){
            double dist2 = DeltaRap2 +std::pow(iphi-phi, 2);
            if (dist2 < d2){
                double cphi = std::cos(iphi), sphi = std::sin(iphi);
                double P0 = ( M_PI/2*cs*s.p.t()
                            + cphi*4./3.*s.p.x()
                            + sphi*4./3.*s.p.y()
                            )*dOmega;              
                pT += 3*P0*(cphi*cosjet+sphi*sinjet);
            }
        }
    }
 }
    for (int j=0; j<results.size(); j++){
        auto jetP = results[j];                  
        f << corrected_pT[j] << " " << jetP.xT() << " " << jetP.phi() << " " 
          << jetP.pseudorap() << " " << jetP.m2() << " " 
          << sigma_gen << std::endl;
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
        double cosjet = std::cos(phi);
        double sinjet = std::sin(phi);

        // The source list is not part of the jet finding, but used as a 
        // correlated back-ground (may be a confusing name) 
        // to further correct the jet energy.
        double cs = 1/1.7;
        int Nphi = 10;
        double dphi = 2.*jetRadius/Nphi;
        double dOmega = 1./4./M_PI*dphi;
        for (auto & s : jlist) {
            // correct jet energy if the source 
            // current point inside the jet radius
            double DeltaRap2 = std::pow(rap - s.x.rap(),2);
            if (DeltaRap2 >= d2) continue;
            for (double iphi=phi-jetRadius; 
                iphi<=phi+jetRadius; iphi+=dphi){
                double dist2 = DeltaRap2 +std::pow(iphi-phi, 2);
                if (dist2 < d2){
                    double cphi = std::cos(iphi), sphi = std::sin(iphi);
	            double P0 = ( M_PI/2*cs*s.p.t()
	                        + cphi*4./3.*s.p.x()
	                        + sphi*4./3.*s.p.y()
                                )*dOmega;              
                    pT += 3*P0*(cphi*cosjet+sphi*sinjet);
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
                double r = std::sqrt(std::pow(phi-p.p.phi(),2)
                                    +std::pow(rap-p.p.pseudorap(),2));
                if (r < 3){
                    int index = 0;
                    for (index=0; index<17; index++){
                        if (r>=bins[index] && r<bins[index+1]) break;
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
            double cs = 1/1.7;
            int Nphi = 300;
            double dphi = 2.*3/Nphi;
            double dOmega = 1./4./M_PI*dphi;
            for (auto & s : jlist) {
                // correct jet energy if the source 
                // current point inside the jet radius
                double DeltaRap2 = std::pow(rap - s.x.rap(),2);
                if (DeltaRap2 >= 3*3) continue;
                for (double iphi=phi-3; 
                    iphi<=phi+3; iphi+=dphi){
                    double dist = std::sqrt(DeltaRap2+std::pow(iphi-phi, 2));
                    if (dist < 3.){
                        double cphi = std::cos(iphi), sphi = std::sin(iphi);
	                double P0 = ( M_PI/2*cs*s.p.t()
	                            + cphi*4./3.*s.p.x()
	                            + sphi*4./3.*s.p.y()
                                    )*dOmega;              
                        int index = 0;
                        for (index=0; index<17; index++){
                            if (dist>=bins[index] && dist<bins[index+1]) break;
                        }
                        double jpT = 3*P0;
                        hist[index] += jpT;
                        w1 += jpT;
                    }   
                }
            }

            for (int i=0; i<hist.size(); i++) 
                f << hist[i]/(bins[i+1]-bins[i]) << " "; 
            f << w1 << " ";
            for (int i=0; i<hist0.size(); i++)
                f << hist0[i]/(bins[i+1]-bins[i]) << " "; 
            f << w0 << " ";
            f << sigma_gen << std::endl;
        }
    }
}


void JetShapeTower(std::vector<particle> plist, 
             std::vector<current> jlist,
             double jetRadius, 
             double jetpTMin, 
             double jetyMin, 
             double jetyMax,
             std::string fname, double sigma_gen){
    // contruct four momentum tower in the eta-phi plane
    int Neta = 200, Nphi = 200;
    double etamin = -3., etamax = 3.;
    double deta = (etamax-etamin)/Neta; 
    // Phi range has to be the same as the range of std::atan2(*,*);
    double phimin = -M_PI, phimax = M_PI; 
    double dphi = (phimax-phimin)/Nphi;
    double dOmega = 1./4./M_PI*dphi;
    std::vector<std::vector<fourvec> > towers;
    towers.resize(Neta);
    for (auto & it : towers) {
        it.resize(Nphi);
        for (auto & iit: it) iit = fourvec{0.,0.,0.,0.};
    }
    // put hard particles into the tower
    for (auto & p : plist){
        double eta = p.p.pseudorap();
        double phi = p.p.phi();
        int ieta = corp_index(eta, etamin, etamax, deta, Neta);
        if (ieta<0) continue;
        int iphi = corp_index(phi, phimin, phimax, dphi, Nphi);
        towers[ieta][iphi] = towers[ieta][iphi] + p.p;
    }
    // source integral contribution to each tower

    for (auto & s : jlist) {
        int ieta = corp_index(s.x.rap(), etamin, etamax, deta, Neta);
        if (ieta<0) continue;
        for (int iphi=0; iphi<Nphi; iphi++){
            double phi = phimin+iphi*dphi;
            double cphi = std::cos(phi), sphi = std::sin(phi);
            double P0 = ( M_PI/2.*s.cs*s.p.t()
                        + cphi*4./3.*s.p.x()
                        + sphi*4./3.*s.p.y()
                        )*dOmega;              
            fourvec k{1., s.cs*cphi, s.cs*sphi, 0};   
            double newphi = k.boost_back(s.v[0], s.v[1], 0).phi();   
            int newiphi = corp_index(newphi, phimin, phimax, dphi, Nphi);
            fourvec dpmu{P0/s.cs, 3*P0*cphi, 3*P0*sphi, 0};
            towers[ieta][newiphi] = towers[ieta][newiphi]+
                                dpmu.boost_back(s.v[0], s.v[1], 0).boost_back(0,0, s.v[2]);
        }
    }
    // Use the towers to do jet finding: anti-kT
    int power = -1; 
    // Define a jet
    //for (int i=0; i<Rs.size(); i++){

    //double jetRadius = Rs[i];
    fastjet::JetDefinition jetDef(fastjet::genkt_algorithm, jetRadius, power);
    // Objetcts to operate on
    std::vector<fastjet::PseudoJet> fjInputs;
    // cuts
    fastjet::Selector select_rapidity = 
          fastjet::SelectorRapRange(jetyMin, jetyMax);
    // when we do the first jet finding: 
    // set those towers with negative contribution to zero.
    for (int ieta=0; ieta<Neta; ieta++) {
        double eta = etamin+ieta*deta;
        double kz = std::sinh(eta), k0 = std::cosh(eta);
        for (int iphi=0; iphi<Nphi; iphi++) {
            auto iit = towers[ieta][iphi];
            double phi = phimin+iphi*dphi;
            double kx = std::cos(phi), ky = std::sin(phi);
            double pT = iit.x()*kx + iit.y()*ky;
            if (iit.t()>0 && pT>0){
                fastjet::PseudoJet fp(iit.x(), iit.y(), iit.z(), iit.t());
                fp.set_user_info(new MyInfo(0, 0, sigma_gen));
                fjInputs.push_back(fp);
            }
        }
    }
    std::vector<fastjet::PseudoJet> inclusiveJets, jets;
    fastjet::ClusterSequence clustSeq(fjInputs, jetDef);
    inclusiveJets = clustSeq.inclusive_jets(jetpTMin);
    std::vector<fastjet::PseudoJet> selected_jets = select_rapidity(inclusiveJets);
    jets = sorted_by_pt(selected_jets);
     
    std::ofstream f(fname.c_str(), std::ios_base::app);
    // once the jet is find, recombine all the bins with negative contribution
    std::vector<fourvec> results;
    std::vector<double> pT0;
    for (auto & j : jets){
        fourvec jetP{j.e(), j.px(), j.py(), j.pz()};
        fourvec newP{0., 0., 0., 0.};
        pT0.push_back(jetP.xT());
        double Jphi = jetP.phi();
        double Jeta = jetP.pseudorap();
        for (double eta=Jeta-jetRadius; eta<=Jeta+jetRadius; eta+=deta){
            int ieta = corp_index(eta, etamin, etamax, deta, Neta);
            if (ieta<0) continue;
            double Rp = std::sqrt(1e-9+std::pow(jetRadius,2)-std::pow(Jeta-eta,2)); 

            for (double phi=Jphi-Rp; phi<=Jphi+Rp; phi+=dphi){
                double newphi = phi;
                if (phi<-M_PI) newphi+=2*M_PI;
                if (phi>M_PI) newphi-=2*M_PI;
                int iphi = corp_index(newphi, phimin, phimax, dphi, Nphi);
                newP = newP + towers[ieta][iphi];
            }  
        }
        results.push_back(newP);
        double bins[18] = {0, .05, .1, .15, .2, .25, .3,
                        .35, .4, .45, .5, .6, .7, .8, 1,
                        1.5, 2.0, 3.0};
        if (newP.xT()>120 && std::abs(newP.pseudorap())<1.6){
            std::vector<double> hist;
            hist.resize(17);
            for (auto & it : hist) it = 0.;

            for (double eta=Jeta-2.; eta<=Jeta+2.; eta+=deta){

                int ieta = corp_index(eta, etamin, etamax, deta, Neta);
                if (ieta<0) continue;
                double Rp = std::sqrt(1e-9+2*2
                                          -std::pow(Jeta-eta,2)); 
                for (double phi=Jphi-Rp; phi<=Jphi+Rp; phi+=dphi){
                    double dist = std::sqrt(std::pow(Jeta-eta,2)
                                          + std::pow(Jphi-phi,2));
                    double newphi = phi;
                    int index = 0;
                    for (index=0; index<17; index++){
                        if (dist>=bins[index] && dist<bins[index+1]) break;
                    }
                    if (phi<-M_PI) newphi+=2*M_PI;
                    if (phi>M_PI) newphi-=2*M_PI;
                    int iphi = corp_index(newphi, phimin, phimax, dphi, Nphi);
                    hist[index] += towers[ieta][iphi].x()*std::cos(phi)+towers[ieta][iphi].y()*std::sin(phi);
                }  
            }
            for (int i=0; i<hist.size(); i++) 
                f << hist[i]/(bins[i+1]-bins[i]) << " "; 
            f << sigma_gen << std::endl;
        }
    }    
}

double F(int order, double vperp, double dphi, double deta){
    double g = 1./std::sqrt(1-vperp*vperp);
    
    return std::pow(g*(std::cosh(deta)-vperp*std::cos(dphi)), order);
}

double Fa(double a){
    double d = std::sqrt(1.-a*a);
    return 4*std::atan((a+1)/d)/d;
}

double Norm(int order, double vperp, double dphi){
    double x = vperp*std::cos(dphi);
    double dx = std::max(0.1*x, 0.02);
    if (order==3){
        return (1./8*Fa(x-3*dx) - Fa(x-2*dx) + 13./8.*Fa(x-dx)
              - 1./8*Fa(x+3*dx) + Fa(x+2*dx) - 13./8.*Fa(x+dx)
               )/std::pow(dx, 3);
    }
    if (order==4){
        return (- 1./6*Fa(x-3*dx) + 2.*Fa(x-2*dx) - 13./2.*Fa(x-dx)
                + 28./3.*Fa(x)
                + 1./6*Fa(x+3*dx) - 2.*Fa(x+2*dx) + 13./2.*Fa(x+dx)
               )/std::pow(dx, 4);
    }
    else{
        LOG_INFO << "should not do this!";
        return 0.;
    }
}
void TestSource(
             std::vector<current> jlist,
             std::string fname) {
    // contruct four momentum tower in the eta-phi plane
    int Neta = 43, Nphi = 48;
    double etamin = -3.15, etamax = 3.15;
    double deta = (etamax-etamin)/Neta; 
    // Phi range has to be the same as the range of std::atan2(*,*);
    double phimin = -M_PI, phimax = M_PI; 
    double dphi = (phimax-phimin)/Nphi;
    std::vector<std::vector<fourvec> > towers;
    towers.resize(Neta);
    std::vector<std::vector<double> > dpTarray;
    towers.resize(Neta);
    dpTarray.resize(Neta);
    for (auto & it : towers) {
        it.resize(Nphi);
        for (auto & iit: it) iit = fourvec{0.,0.,0.,0.};
    }
    for (auto & it : dpTarray) {
        it.resize(Nphi);
        for (auto & iit: it) iit = 0;
    }
    // source integral contribution to each tower
    for (auto & s : jlist) {
        double cs = s.cs;
        double srap = s.x.rap();
        for (int iphi=0; iphi<Nphi; iphi++){
            double phi0 = phimin+iphi*dphi;
            double cphi0 = std::cos(phi0), sphi0 = std::sin(phi0);
              
            for (double w=-3; w<=3; w+=.2){    
              double costheta = std::tanh(w);  
              double sintheta = std::sqrt(1.-costheta*costheta);  
              double P0 = (cs*s.p.t()
                        + cphi0*sintheta*s.p.x()
                        + sphi0*sintheta*s.p.y()
                        + costheta*s.p.z()
                        )*1./4./M_PI*dphi*.2/std::pow(std::cosh(w),2);
              fourvec dpmu{P0/cs,
                           3*P0*cphi0*sintheta, 
                           3*P0*sphi0*sintheta, 
                           3*P0*costheta};

              fourvec k{1, cs*sintheta*cphi0, cs*sintheta*sphi0, cs*costheta};              
              k = k.boost_back(s.v[0], s.v[1], 0.0);
              double phi_v0 = std::atan2(s.v[1], s.v[0]);
              
              double vcx = 0.7*std::cos(k.phi());
              double vcy = 0.7*std::sin(k.phi());
              double vperp = std::sqrt(vcx*vcx+vcy*vcy);
              double gamma_perp = 1./std::sqrt(1-vperp*vperp);

              double finaletas = s.x.rap()+k.rap();

              fourvec Umu={gamma_perp*std::cosh(finaletas), 
                          gamma_perp*vcx,
                          gamma_perp*vcy,  
                          gamma_perp*std::sinh(finaletas)};
              double phiv = Umu.phi();
              dpmu = dpmu.boost_back(s.v[0], s.v[1], 0).boost_back(0, 0, s.v[2]);

              double UdotdP = dot(Umu, dpmu);
              for (int i=0; i<Nphi*2.; i++){
                  double deltaphi = .1+phimin+i*dphi/2.;
                  double finalphi = phiv+deltaphi;
                  if (finalphi<-M_PI) finalphi+=2*M_PI;
                  if (finalphi>M_PI) finalphi-=2*M_PI;

                  for (double dw=-6; dw<=6; dw+=.15){
                   double newfinaletas = finaletas + dw;
                   int ieta = corp_index(newfinaletas, etamin, etamax, deta, Neta);
                   if (ieta<0) continue;

                   fourvec Nmu{std::cosh(newfinaletas), 
                               std::cos(finalphi), 
                               std::sin(finalphi), 
                               std::sinh(newfinaletas)};
                   int finalphiindex = corp_index(finalphi, phimin, phimax, dphi, Nphi);
                   double prefactor = std::pow(M_PI, 2)/4./std::pow(2*M_PI, 3);
                   double A = prefactor*42/F(3, vperp, deltaphi, dw);
                   double B = -prefactor*18/F(4, vperp, deltaphi, dw);
                   double deltapT = (A*UdotdP + B*dot(Nmu, dpmu))*dphi/2.*.15;
                   dpTarray[ieta][finalphiindex] += deltapT;
                  }
              }
            }
        }
    }
    std::ofstream f(fname.c_str());
    for (int ieta=0; ieta<Neta; ieta++) {
        double eta = etamin+ieta*deta;
        for (int iphi=0; iphi<Nphi; iphi++){
            double phi = phimin+iphi*dphi;
            double cphi = std::cos(phi), sphi = std::sin(phi);
            auto it = dpTarray[ieta][iphi];
            f << it << " ";
        }
        f << std::endl;
    }   
}

