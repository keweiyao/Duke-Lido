#include <fstream>
#include <string>
#include <math.h>
#include <vector>
#include "jet_finding.h"
#include "workflow.h"
#include "lorentz.h"
#include "simpleLogger.h"
#include "integrator.h"
#include <sstream>

inline int corp_index(double x, double xL, double xH, double dx, int Nx){
    if (x<xL || x>xH) return -1;
    else if (x==xH) return Nx-1;
    else return int((x-xL)/dx);
}

double inline F(double gamma_perp, double vperp, double dphi, double deta){
    return gamma_perp*(std::cosh(deta)-vperp*std::cos(dphi));
}

void redistribute(
          std::vector<current> & jlist, 
          std::vector<std::vector<fourvec> > & DeltaPmu,
          std::vector<std::vector<double> > & DeltaPT,
          const int _Ny, const int _Nphi, 
          const double ymin, const double ymax,
          const double phimin, const double phimax,
          int coarse_level){
    // use coarse grid
    int Ny = int(_Ny/coarse_level);
    int Nphi = int(_Nphi/coarse_level);
    const double dy = (ymax-ymin)/(Ny-1); 
    const double dphi = (phimax-phimin)/Nphi;
    const double prefactor = std::pow(M_PI, 2)/4./std::pow(2*M_PI, 3)/4./M_PI;
    const double vradial = 0.6;
    const double gamma_radial = 1./std::sqrt(1-std::pow(vradial,2));

    std::vector<std::vector<fourvec> > CoarsedPmu;
    std::vector<std::vector<double> > CoarsedPT;
    CoarsedPmu.resize(Ny);
    for (auto & it : CoarsedPmu) {
        it.resize(Nphi);
        for (auto & iit: it) iit = fourvec{0.,0.,0.,0.};
    }
    CoarsedPT.resize(Ny);
    for (auto & it : CoarsedPT) {
        it.resize(Nphi);
        for (auto & iit: it) iit = 0.;
    }


    LOG_INFO << jlist.size() << " sources";
    for (int iy=0; iy<Ny; iy++){
        LOG_INFO << iy << " " << ymin+iy*dy;
        double y = ymin+iy*dy;    
        double chy = std::cosh(y), 
               shy = std::sinh(y);
        for (int iphi=0; iphi<Nphi; iphi++){
            double phi = phimin+(iphi+.5)*dphi;
            double cphi = std::cos(phi), sphi = std::sin(phi);
            fourvec Nmu0{chy,cphi,sphi,shy};
            auto code = [jlist,
                         phi, cphi, sphi, 
                         y, chy, shy, 
                         vradial,gamma_radial]
                         (const double * X){
                std::vector<double> res;
                res.push_back(0.);
                double cthetak = X[0], phik = X[1];
                double sthetak = std::sqrt(1.-std::pow(cthetak,2));
                double cphik=std::cos(phik), 
                       sphik=std::sin(phik);
                for (auto & source : jlist){
                    double etas = source.rap;
                    fourvec k{1./source.cs, 
                              sthetak*cphik,
                              sthetak*sphik, 
                              cthetak};     
                    double vradial_x = vradial*k.x()/k.xT();
                    double vradial_y = vradial*k.y()/k.xT();
                    double krap = k.rap();                   
                    double G0 = source.p.t()*source.cs
                              + source.p.x()*k.x()
                              + source.p.y()*k.y()
                              + source.p.z()*k.z();
                    double G03 = 3.*G0;
                    fourvec dgmu = {G0/source.cs,
                                    G03*k.x(), 
                                    G03*k.y(), 
                                    G03*k.z()};
                    double kT = k.xT();

                    double chketa = (k.t()*source.x.t()+
                              k.z()*source.x.z())/source.tau/k.tau();
                    double shketa = (k.z()*source.x.t()+
                              k.t()*source.x.z())/source.tau/k.tau();
                    double UdotdG = gamma_radial*(
                             chketa*dgmu.t() - shketa*dgmu.z() 
                         - vradial_x*dgmu.x() - vradial_y*dgmu.y()
                               );
                    double NdotdG = chy*dgmu.t() - shy*dgmu.z()
                                  - cphi*dgmu.x() - sphi*dgmu.y();
                    double Ndotdk = chy*k.t() - shy*k.z()
                                  - cphi*k.x() - sphi*k.y();
                    double NdotU = gamma_radial*(
                             chketa*chy - shketa*shy
                         - vradial_x*cphi - vradial_y*sphi
                               );
                    double sigma = gamma_radial*(
                                    chy*chketa - shy*shketa
                               - cphi*vradial_x - sphi*vradial_y
                                   );
                    res[0] += (32.*UdotdG*sigma-24.*NdotdG)/std::pow(sigma,4);
                }
                return res;
            };
	    double xmin[2] = {-1., -M_PI};
	    double xmax[2] = {1., M_PI};
            double error;
            double res = quad_nd(code, 2, 1, xmin, xmax, 
                                error, 0., .05, 100)[0] 
                       * prefactor;
            CoarsedPmu[iy][iphi] = CoarsedPmu[iy][iphi] + Nmu0*res;
            CoarsedPT[iy][iphi] += res;
        }
    }
    // Interpolate the coarse grid to finer grid
    const double fine_dy = (ymax-ymin)/(_Ny-1); 
    const double fine_dphi = (phimax-phimin)/_Nphi;
    for (int i=0; i<_Ny; i++){
        double y = ymin + fine_dy*i;
        int ii = corp_index(y, ymin, ymax, dy, Ny);
        if (ii<0) continue;
        if (ii==Ny-1) ii--;
        double residue1 = (y-(ymin+dy*ii))/dy;
        double u[2] = {1.-residue1, residue1};
        for (int j=0; j<_Nphi; j++){ 
            double phi = phimin + fine_dphi*(j+.5);
            int jj = corp_index(phi, phimin, phimax, dphi, Nphi);
            if (jj==Nphi-1) jj--;
            double residue2 = (phi-(phimin+dphi*(jj+.5)))/dphi;
            double v[2] = {1.-residue2, residue2};
            for (int k1=0; k1<2; k1++){
                for (int k2=0; k2<2; k2++){
                    LOG_INFO << "a"<< DeltaPmu[i][j] << " " << CoarsedPmu[ii+k1][jj+k2];
                    DeltaPmu[i][j] = DeltaPmu[i][j] + 
                        CoarsedPmu[ii+k1][jj+k2]*fine_dy*fine_dphi*u[k1]*v[k2];
                    DeltaPT[i][j] += 
                        u[k1]*v[k2]*CoarsedPT[ii+k1][jj+k2]*fine_dy*fine_dphi;
                    LOG_INFO << "b"<< DeltaPmu[i][j];
                }
            }
        }
    }
}

void FindJetTower(std::vector<particle> plist, 
             std::vector<current> jlist,
             std::vector<double> Rs,
             double jetpTMin, 
             double jetyMin, 
             double jetyMax,
             std::string fheader, 
             double sigma_gen) {
    // contruct four momentum tower in the eta-phi plane
    int Neta = 240, Nphi = 240;
    double etamin = -3., etamax = 3.;
    double deta = (etamax-etamin)/(Neta-1); 
    // Phi range has to be the same as the range of std::atan2(*,*);
    double phimin = -M_PI, phimax = M_PI; 
    double dphi = (phimax-phimin)/Nphi;
    std::vector<std::vector<fourvec> > Pmutowers;
    std::vector<std::vector<double> > PTtowers;
    Pmutowers.resize(Neta);
    for (auto & it : Pmutowers) {
        it.resize(Nphi);
        for (auto & iit: it) iit = fourvec{0.,0.,0.,0.};
    }
    PTtowers.resize(Neta);
    for (auto & it : PTtowers) {
        it.resize(Nphi);
        for (auto & iit: it) iit = 0.;
    }

    // put hard particles into the towers
    for (auto & p : plist){
        double eta = p.p.pseudorap();
        double phi = p.p.phi();
        int ieta = corp_index(eta, etamin, etamax, deta, Neta);
        if (ieta<0) continue;
        int iphi = corp_index(phi, phimin, phimax, dphi, Nphi);
        Pmutowers[ieta][iphi] = Pmutowers[ieta][iphi] + p.p;
        PTtowers[ieta][iphi] += p.p.xT();
    }
    // put soft energy-momentum deposition into the towers
    if (jlist.size() >0){
    redistribute(
          jlist, 
          Pmutowers,
          PTtowers,
          Neta, Nphi,
          etamin, etamax,
          phimin, phimax,
          24);
    }
    // Use the towers to do jet finding: anti-kT
    int power = -1; 
    // loop over jetradius
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
                auto ipmu = Pmutowers[ieta][iphi];
                auto ipT = PTtowers[ieta][iphi];
                // when we do the first jet finding,
                // only use towers with positive contribution
                if (ipmu.t()>0 && ipT>0){
                    fastjet::PseudoJet fp(ipmu.x(), ipmu.y(),
                                          ipmu.z(), ipmu.t());
                    fp.set_user_info(new MyInfo(0, 0, 0));
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
        std::stringstream filename_jets;
        filename_jets << fheader << "-jets.dat";
        std::ofstream f1(filename_jets.str(), std::ios_base::app);
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
                                          -std::pow(Jeta-eta,2) ); 
                for (double phi=Jphi-Rp; phi<=Jphi+Rp; phi+=dphi){
                    double newphi = phi;
                    if (phi<-M_PI) newphi+=2*M_PI;
                    if (phi>M_PI) newphi-=2*M_PI;
                    int iphi = corp_index(newphi, phimin, phimax, dphi, Nphi);
                    newP = newP + Pmutowers[ieta][iphi];
                }  
            }
            results.push_back(newP);    
            // label the jet flavor:
            // if heavy, label it by heavy quark id
            // else, label it by the hardest (pT) parton inside
            int flavor=0;
            fourvec pf{0,0,0,0};
            for (auto & p : plist){
                if ( std::abs(p.pid) == 5){
                    double eta = p.p.pseudorap();
                    double phi = p.p.phi();
                    if (std::sqrt(std::pow(eta-Jeta,2)
                                 +std::pow(phi-Jphi,2)) < jetRadius) {
                        flavor = std::abs(p.pid);
                        pf = p.p;
                        break;
                    }
                }
            }      
            if (flavor==0)  {
                for (auto & p : plist){
                    double eta = p.p.pseudorap();
                    double phi = p.p.phi();
                    if (std::sqrt(std::pow(eta-Jeta,2)
                                 +std::pow(phi-Jphi,2)) < jetRadius) {
                        if ( p.p.xT()>pf.xT() ){
                            flavor = std::abs(p.pid);
                            pf = p.p;
                        }
                    }
                }
            }
            f1 << iR << " " 
              << newP.xT() << " " << newP.phi() << " " 
              << newP.pseudorap() << " " << newP.m2() << " " 
              << sigma_gen << " " 
              << flavor << " " 
              << pf.xT() << " " << pf.phi() << " " 
              << pf.pseudorap() << " "  << std::endl;
        }

        // also compute jetshape, with the following bins
        double bins[19] = {0, .05, .1, .15,  .2, .25, .3,
                          .35, .4, .45, .5,  .6, .7,  .8, 
                           1., 1.5, 2.0, 2.5, 3.0};
        std::stringstream filename_jetshape;
        filename_jetshape << fheader << "-jetshape.dat";
        std::ofstream f2(filename_jetshape.str(), std::ios_base::app);
        for (auto newP : results){
            double Jphi = newP.phi();
            double Jeta = newP.pseudorap();
            std::vector<double> hist;
            hist.resize(18);
            for (auto & it : hist) it = 0.;
            for (double eta=Jeta-3.; eta<=Jeta+3.; eta+=deta){

                int ieta = corp_index(eta, etamin, etamax, deta, Neta);
                if (ieta<0) continue;
                double Rp = std::sqrt(1e-9+std::pow(3.,2)
                                          -std::pow(Jeta-eta,2) ); 
                for (double phi=Jphi-Rp; phi<=Jphi+Rp; phi+=dphi){
                    double dist = std::sqrt(std::pow(Jeta-eta,2)
                                          + std::pow(Jphi-phi,2));
                    int index = 0;
                    for (index=0; index<18; index++){
                        if (dist>=bins[index] && dist<bins[index+1]) break;
                    }
                    double newphi = phi;
                    if (phi<-M_PI) newphi+=2*M_PI;
                    if (phi>M_PI) newphi-=2*M_PI;
                    int iphi = corp_index(newphi, phimin, phimax, dphi, Nphi);
                    hist[index] += PTtowers[ieta][iphi];
                }  
            }
            f2 << iR << " " << newP.xT() << " " << newP.phi() << " " 
               << newP.pseudorap() << " " << newP.m2() << " " 
               << sigma_gen << " ";
            for (int i=0; i<hist.size(); i++) 
                f2 << hist[i]/(bins[i+1]-bins[i]) << " "; 
            f2 << std::endl;
        }
    }   
}

void TestSource(
             std::vector<current> jlist,
             std::string fname) {
    // contruct four momentum tower in the eta-phi plane
    int Neta = 43, Nphi = 48;
    double etamin = -3.15, etamax = 3.15;
    double deta = (etamax-etamin)/(Neta-1); 
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

    fourvec T0{0,0,0,0};
    for (auto& r : towers){
        for (auto& it: r){
            T0 = T0 + it;
        }
    }
    if (jlist.size() >0){
    redistribute(
          jlist, 
          towers,
          dpTarray,
          Neta, Nphi,
          etamin, etamax,
          phimin, phimax, 1);
    }

    fourvec T1{0,0,0,0};
    for (auto& r : towers){
        for (auto& it: r){
            T1 = T1 + it;
        }
    }
    LOG_INFO << "total" << T1-T0;

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

