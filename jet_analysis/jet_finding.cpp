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
#include <thread>
#include <gsl/gsl_sf_gamma.h>

bool compare_jet(Fjet A, Fjet B){
    return (A.pT > B.pT);
}

MediumResponse::MediumResponse(std::string _name): name(_name){
    std::vector<size_t> shape;
    std::vector<double> low, high;
    shape.resize(4);
    low.resize(4);
    high.resize(4);
    low[0] = -3.; low[1] = -M_PI; low[2] = 0.; low[3] = 0.0;
    high[0] = 3.; high[1] = M_PI; high[2] = 30.; high[3] = 0.9;
    shape[0] = 30; shape[1] = 30; shape[2] = 31; shape[3] = 10;
    Gmu = std::make_shared<TableBase<fourvec, 4>>(name, shape, low, high);
}

void MediumResponse::load(std::string fname){
    Gmu->Load(fname);
}

void MediumResponse::init(std::string fname){
    LOG_INFO << fname << " Generating tables for approx medium response functions";
    auto code = [this](int start, int end) { this->compute(start, end); };
    std::vector<std::thread> threads;
    size_t nthreads = std::thread::hardware_concurrency();
    size_t padding = size_t(std::ceil(Gmu->length()*1./nthreads));
    for(auto i=0; i<nthreads; ++i) {
        int start = i*padding;
        int end = std::min(padding*(i+1), Gmu->length());
        threads.push_back( std::thread(code, start, end) );
    }
    for(auto& t : threads) t.join();
    Gmu->Save(fname);
}

void MediumResponse::compute(int start, int end){
    std::vector<size_t> index;
    index.resize(4);
    for(auto i=start; i<end; ++i){
        size_t q = i;
        for(int d=4-1; d>=0; d--){
            size_t dim = Gmu->shape(d);
            size_t n = q%dim;
            q = q/dim;
            index[d] = n;
        }
        auto params = Gmu->parameters(index);
        double rap = params[0];
        double phi = params[1];
        double pTmin_over_T = params[2];
        double vperp = params[3];
        double cs = std::sqrt(0.2);
        double chrap = std::cosh(rap);
        double shrap = std::sinh(rap);
        auto code = [rap, chrap, shrap, phi, pTmin_over_T, vperp, cs]
                 (const double * X){
            double costhetak = X[0];
            double phik = X[1];
            double yk = 0.5*std::log((1./cs+costhetak)/(1./cs-costhetak));
            double sinthetak = std::sqrt(1.-costhetak*costhetak);
            double gamma_vperp = 1./std::sqrt(1.-vperp*vperp);
            double chyk = std::cosh(yk), shyk = std::sinh(yk);
            double sigma = gamma_vperp * ( std::cosh(rap-yk)
                                 - vperp * std::cos(phi-phik) );
            double GInc = gsl_sf_gamma_inc_Q(5., pTmin_over_T*sigma);
            double sigma_3rd = std::pow(sigma, 3);
            double sigma_4th = std::pow(sigma, 4);
            double A = (
                     4./3.*gamma_vperp/sigma_3rd * (
                chyk/cs - 3*(sinthetak*vperp + costhetak*shyk)
                     )
                   - 1./sigma_4th * (
                chrap/cs - 3*(sinthetak*std::cos(phi-phik)+ shrap*costhetak)
                     )
                  ) * GInc * 3. / std::pow(4.*M_PI, 2);
            std::vector<double> res;
            res.resize(4);
            res[0] = A*cs;
            res[1] = A*sinthetak*std::cos(phik);
            res[2] = A*sinthetak*std::sin(phik);
            res[3] = A*costhetak;
            return res;
        };
        double wmin[2] = {-.999, -M_PI};
        double wmax[2] = {.999, M_PI};
        double error;
        std::vector<double> res = quad_nd(code, 2, 4, wmin, wmax, error, 1e-7);
        fourvec gmu{res[0], res[1], res[2], res[3]};
        Gmu->SetTableValue(index, gmu);
    }
}

double MediumResponse::get_dpT_dydphi(double rap, double phi, fourvec Pmu, double vperp, double pTmin_over_T){
    auto gmu = Gmu->InterpolateTable({rap, phi, pTmin_over_T, vperp});
    return gmu.t()*Pmu.t()+gmu.x()*Pmu.x()+gmu.y()*Pmu.y()+gmu.z()*Pmu.z();
}

inline int corp_index(double x, double xL, double xH, double dx, int Nx){
    if (x<xL || x>xH) return -1;
    else if (x==xH) return Nx-1;
    else return int((x-xL)/dx);
}

std::vector<Fjet> FindJetTower(
             MediumResponse MR,
             std::vector<particle> plist, 
             std::vector<current> SourceList,
             std::vector<HadronizeCurrent> HadronizeList,
             std::vector<double> Rs, std::vector<double> rbins,
             double jetpTMin, 
             double jetyMin, 
             double jetyMax,
             double sigma_gen, 
	     double pTmin) {
    LOG_INFO << "do jet finding";
    int coarse_level = 5;
    std::vector<Fjet> Results;
    // contruct four momentum tower in the eta-phi plane
    int Neta = 300, Nphi = 300;
    int coarseNeta = int(Neta/coarse_level), 
        coarseNphi = int(Nphi/coarse_level);
    double etamin = -3., etamax = 3.;
    double deta = (etamax-etamin)/(Neta-1); 
    double coarsedeta = (etamax-etamin)/(coarseNeta-1); 
    // Phi range has to be the same as the range of std::atan2(*,*);
    double phimin = -M_PI, phimax = M_PI; 
    double dphi = (phimax-phimin)/Nphi;
    double coarsedphi = (phimax-phimin)/coarseNphi;
    std::vector<std::vector<fourvec> > Pmutowers;
    Pmutowers.resize(Neta);
    for (auto & it : Pmutowers) {
        it.resize(Nphi);
        for (auto & iit: it) iit = fourvec{0.,0.,0.,0.};
    }
    std::vector<std::vector<double> > PTtowers;
    PTtowers.resize(Neta); 
    for (auto & it : PTtowers) {
        it.resize(Nphi);
        for (auto & iit: it) iit = 0.;
    }
    std::vector<std::vector<std::vector<double> > > coarsePT;
    coarsePT.resize(2);
    for (auto & t : coarsePT) {
       t.resize(coarseNeta); 
       for (auto & it : t) {
            it.resize(coarseNphi);
            for (auto & iit: it) iit = 0.;
       }
    }

    // put hard particles into the towers
    for (auto & p : plist){
        double eta = p.p.pseudorap();
        double phi = p.p.phi();
        int ieta = corp_index(eta, etamin, etamax, deta, Neta);
        if (ieta<0) continue;
        int iphi = corp_index(phi, phimin, phimax, dphi, Nphi);
        Pmutowers[ieta][iphi] = Pmutowers[ieta][iphi] + p.p;
        if (p.p.xT()>pTmin)
            PTtowers[ieta][iphi] += p.p.xT();
    }
    // put soft energy-momentum deposition into the towers
    std::vector<current> clist;
    clist.resize(coarseNeta);
    for (int ieta=0; ieta<coarseNeta; ieta++){
        clist[ieta].chetas = std::cosh(etamin+ieta*coarsedeta);
        clist[ieta].shetas = std::sinh(etamin+ieta*coarsedeta);
        clist[ieta].p.a[0] = 0.;
        clist[ieta].p.a[1] = 0.;
        clist[ieta].p.a[2] = 0.;
        clist[ieta].p.a[3] = 0.;
    }
    for (auto & s: SourceList){
        int i = corp_index(std::atanh(s.shetas/s.chetas), etamin, etamax, coarsedeta, coarseNeta);
        if (i<0) continue;
        clist[i].p = clist[i].p + s.p;
    }
    for (int ieta=0; ieta<coarseNeta; ieta++) {
        double eta = etamin+ieta*coarsedeta;
        for (int iphi=0; iphi<coarseNphi; iphi++) {
            double phi = phimin+iphi*coarsedphi;
            double dpT = 0., dpTcut = 0.;
            for (auto & s: clist){
                double etas = std::atanh(s.shetas/s.chetas);
                dpT += MR.get_dpT_dydphi(eta-etas, phi, s.p, 0.6, 0.);
                dpTcut += MR.get_dpT_dydphi(eta-etas, phi, s.p, 0.6, pTmin/0.165);
            }
            coarsePT[0][ieta][iphi] = dpT;
            coarsePT[1][ieta][iphi] = dpTcut;
        }
    }
    for (int ieta=0; ieta<Neta; ieta++) {
        double eta = etamin+ieta*deta;
        double chy = std::cosh(eta), shy = std::sinh(eta);
        int ii = corp_index(eta, etamin, etamax, coarsedeta, coarseNeta);
        if (ii<0) continue;
        if (ii==coarseNeta-1) ii--;
        double residue1 = (eta-(etamin+coarsedeta*ii))/coarsedeta;
        double u[2] = {1.-residue1, residue1};
        for (int iphi=0; iphi<Nphi; iphi++) {
            double phi = phimin+iphi*dphi;
            double cphi = std::cos(phi), sphi = std::sin(phi);
            fourvec Nmu0{chy, cphi, sphi, shy};
            int jj = corp_index(phi, phimin, phimax, coarsedphi, coarseNphi);
            if (jj==coarseNphi-1) jj--;
            double residue2 = (phi-(phimin+coarsedphi*(jj)))/coarsedphi;
            double v[2] = {1.-residue2, residue2};
            for (int k1=0; k1<2; k1++){
                for (int k2=0; k2<2; k2++){
                    Pmutowers[ieta][iphi] = Pmutowers[ieta][iphi] + 
                        Nmu0*coarsePT[0][ii+k1][jj+k2]*deta*dphi*u[k1]*v[k2];
                    PTtowers[ieta][iphi] += 
                        coarsePT[1][ii+k1][jj+k2]*deta*dphi*u[k1]*v[k2];
                }
            }
        }
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
            for (int iphi=0; iphi<Nphi; iphi++) {
                auto ipmu = Pmutowers[ieta][iphi];
                double phi = phimin+iphi*dphi;
                // when we do the first jet finding,
                // only use towers with positive contribution
                if (ipmu.t()>0.){
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
        
        for (auto & j : jets){
            Fjet J;
	    J.sigma = sigma_gen;
            fourvec jetP{j.e(), j.px(), j.py(), j.pz()};
            fourvec newP{0., 0., 0., 0.};
            double Jphi = jetP.phi();
            double Jeta = jetP.pseudorap();

            for (double eta=Jeta-jetRadius; eta<=Jeta+jetRadius; eta+=deta){
                int ieta = corp_index(eta, etamin, etamax, deta, Neta);
                if (ieta<0) continue;
                double Rp = std::sqrt(1e-9+std::pow(jetRadius,2)
                                          -std::pow(Jeta-eta,2) ); 
                double chy = std::cosh(eta), shy = std::sinh(eta);
                for (double phi=Jphi-Rp; phi<=Jphi+Rp; phi+=dphi){
                    double newphi = phi;
                    if (phi<-M_PI) newphi+=2*M_PI;
                    if (phi>M_PI) newphi-=2*M_PI;
                    double cphi = std::cos(newphi), sphi = std::sin(newphi);
                    int iphi = corp_index(newphi, phimin, phimax, dphi, Nphi);
                    fourvec Nmu0{chy,cphi,sphi,shy};
                    newP = newP + Pmutowers[ieta][iphi];
                }  
            }
        J.pmu = newP;
        J.R = jetRadius;
        J.pT = newP.xT();
        J.phi = newP.phi();
        J.eta = newP.pseudorap();
        J.M2 = newP.m2();
            // label the jet flavor:
            // if heavy, label it by heavy quark id
            // else, label it by the hardest (pT) parton inside
            for (auto & p : plist){
                int absid = std::abs(p.pid);
                bool trigger = absid == 4 || 
                             absid == 411 || absid == 421 ||
                             absid == 413 || absid == 423 ||
                             absid == 5 || 
                             absid == 511 || absid == 521 ||
                             absid == 513 || absid == 523;
                if (trigger){
                    double eta = p.p.pseudorap();
                    double phi = p.p.phi();
                    if (std::sqrt(std::pow(eta-Jeta,2)
                        +std::pow(phi-Jphi,2)) < 0.3) { 
                       // for R jet , find heavy flavor in 0.3 bins
                       particle flavor_tag;
                       flavor_tag.p = p.p;
                       flavor_tag.pid = std::abs(p.pid); 
                       J.Ftags.push_back(flavor_tag);
                    }
                }
            } 
        Results.push_back(J);
        }

        for (auto & J : Results){
            double Jphi = J.phi;
            double Jeta = J.eta;

        // compute dpT/dr, makes sense for both pp and AA
        J.shape.resize(rbins.size()-1);
            for (auto & it : J.shape) it = 0.;
            for (double eta=Jeta-3.; eta<=Jeta+3.; eta+=deta){
                int ieta = corp_index(eta, etamin, etamax, deta, Neta);
                if (ieta<0) continue;
                double Rp = std::sqrt(1e-9+std::pow(3.,2)
                                          -std::pow(Jeta-eta,2) ); 
                for (double phi=Jphi-Rp; phi<=Jphi+Rp; phi+=dphi){
                    double dist = std::sqrt(std::pow(Jeta-eta,2)
                                          + std::pow(Jphi-phi,2));
                    int index = 0;
                    for (index=0; index<J.shape.size(); index++){
                        if (dist>=rbins[index] && dist<rbins[index+1]) break;
                    }
                    double newphi = phi;
                    if (phi<-M_PI) newphi+=2*M_PI;
                    if (phi>M_PI) newphi-=2*M_PI;
                    int iphi = corp_index(newphi, phimin, phimax, dphi, Nphi);
                    J.shape[index] += PTtowers[ieta][iphi];
                }  
            }
            for (int i=0; i<J.shape.size(); i++) 
                J.shape[i] /= (rbins[i+1]-rbins[i]);

        // compute dNch/dr, currently only makes sense to pp
            J.dndr.resize(rbins.size()-1);
            for (auto & it : J.dndr) it = 0.;
        for (auto & p : plist){
        if ((!p.charged) || (p.p.xT()<.7)) continue;
        double dist = std::sqrt(std::pow(J.phi-p.p.phi(), 2)
                            + std::pow(J.eta-p.p.pseudorap(), 2) );
        int index = 0;
                for (index=0; index<J.dndr.size(); index++){
                    if (dist>=rbins[index] && dist<rbins[index+1]) break;
                }
                J.dndr[index] += 1.;
            }
            for (int i=0; i<J.dndr.size(); i++)
                J.dndr[i] /= (rbins[i+1]-rbins[i]);
        }
    }
    return Results; 
}

void TestSource(
             MediumResponse MR,
             std::vector<current> SourceList,
             std::vector<HadronizeCurrent> HadronizeList,
             std::string fname) {
    // contruct four momentum tower in the eta-phi plane
    int Neta = 100, Nphi = 101;
    double etamin = -3, etamax = 3;
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

    for (int ieta=0; ieta<Neta; ieta++) {
        double eta = etamin+ieta*deta;
        double chy = std::cosh(eta), shy = std::sinh(eta);
        for (int iphi=0; iphi<Nphi; iphi++) {
            double phi = phimin+iphi*dphi;
            double cphi = std::cos(phi), sphi = std::sin(phi);
            fourvec Nmu0{chy, cphi, sphi, shy};
            double dpT = 0.;
            for (auto & s: SourceList){
                double etas = std::atanh(s.shetas/s.chetas);
                dpT += MR.get_dpT_dydphi(eta, phi, s.p, 0.0, 0.0)*dphi*deta;
            }
            dpTarray[ieta][iphi] = dpT;
            towers[ieta][iphi] = Nmu0*dpT;
        }
    }


    fourvec T1{0,0,0,0};
    for (auto& r : towers){
        for (auto& it: r){
            T1 = T1 + it;
        }
    }
    LOG_INFO << "total " << T1;

    std::ofstream f(fname.c_str());
    for (int ieta=0; ieta<Neta; ieta++) {
        double eta = etamin+ieta*deta;
        for (int iphi=0; iphi<Nphi; iphi++){
            double phi = phimin+iphi*dphi;
            auto it = dpTarray[ieta][iphi];
            f << it << " ";
        }
        f << std::endl;
    }   
}


LeadingParton::LeadingParton(std::vector<double> _pTbins){
    pTbins = _pTbins;
    NpT = pTbins.size()-1;
    
    nchg.resize(NpT); for (auto & it:nchg) it=0.;
    npi.resize(NpT); for (auto & it:npi) it=0.;
    nD.resize(NpT); for (auto & it:nD) it=0.;
    nB.resize(NpT); for (auto & it:nB) it=0.;
    binwidth.resize(NpT); 
    for (int i=0; i<NpT; i++) 
        binwidth[i]=pTbins[i+1]-pTbins[i];
}
void LeadingParton::add_event(std::vector<particle> plist, double sigma_gen){
    for (auto & p : plist){
        int pid = std::abs(p.pid);
        bool ispi = pid==111 || pid == 211;
        bool ischg = pid==111 || pid == 211 || pid == 311 || pid== 321 || pid==2212;
        bool isD = pid==411 || pid == 421 || pid == 413 || pid == 423 || pid ==4;
        bool isB = pid==511 || pid == 521 || pid == 513 || pid == 523 || pid ==5;
        if (std::abs(p.p.rap()) < 1.) {
            double pT = p.p.xT();
            int ii=0;
            for (int i=0; i<NpT; i++){
                if ( (pTbins[i]<pT) && (pT<pTbins[i+1]) ) {
                    ii = i;
                    break;
                }
            }
            if (ii==NpT) ii=NpT-1;
            if (ispi) npi[ii] += sigma_gen;
            if (ischg) nchg[ii] += sigma_gen;
        if (isD) nD[ii] += sigma_gen;
            if (isB) nB[ii] += sigma_gen;
        }
    }
}
void LeadingParton::write(std::string fheader){
    std::stringstream filename;
    filename << fheader << "-LeadingHadron.dat";
    std::ofstream f(filename.str());

    for (int i=0; i<NpT; i++) {
        f    << (pTbins[i]+pTbins[i+1])/2. << " "
             << pTbins[i] << " "
             << pTbins[i+1] << " "
             << nchg[i]/binwidth[i] << " "
             << npi[i]/binwidth[i] << " " 
             << nD[i]/binwidth[i] << " "
             << nB[i]/binwidth[i] << std::endl;
    }
    f.close();
}


JetStatistics::JetStatistics(std::vector<double> _pTbins, 
		std::vector<double> _Rs,
     std::vector<double> _shape_pTbins, 
     std::vector<double> _shape_rbins){
    AllJets.clear();
    for (int i=0; i<11; i++)
	xJbins.push_back(i*0.06);
    xJ.resize(xJbins.size()-1);
    for (auto & it : xJ) it = 0.;
    pTbins = _pTbins;
    Rs = _Rs;
    shape_pTbins = _shape_pTbins;
    shape_rbins = _shape_rbins;
    NpT = pTbins.size()-1;
    shape_NpT = shape_pTbins.size()-1;
    shape_Nr = shape_rbins.size()-1;
    shape_w.resize(shape_NpT);
    for (auto & it:shape_w) it=0.;
    shapes.resize(shape_NpT); 
    for (auto & it:shapes) {
        it.resize(_shape_rbins.size()-1);
        for (auto & iit : it ) iit=0.;
    }
    dnchdr.resize(shape_NpT);
    for (auto & it:dnchdr) {
            it.resize(_shape_rbins.size()-1);
        for (auto & iit : it ) iit=0.;
    }
    dsigmadpT.resize(Rs.size());
    for (auto & it:dsigmadpT) {
            it.resize(NpT);
        for (auto & iit : it ) iit=0.;
    }
    dDdpT.resize(Rs.size());
    for (auto & it:dDdpT) {
            it.resize(NpT);
        for (auto & iit : it ) iit=0.;
    }
    dBdpT.resize(Rs.size());
    for (auto & it:dBdpT) {
            it.resize(NpT);
        for (auto & iit : it ) iit=0.;
    }
    dD0dpT.resize(Rs.size());
    for (auto & it:dD0dpT) {
            it.resize(NpT);
        for (auto & iit : it ) iit=0.;
    }
    dB0dpT.resize(Rs.size());
    for (auto & it:dB0dpT) {
            it.resize(NpT);
        for (auto & iit : it ) iit=0.;
    }


    binwidth.resize(NpT);
    for (int i=0; i<NpT; i++)
        binwidth[i]=pTbins[i+1]-pTbins[i];
}
void JetStatistics::add_event(std::vector<Fjet> jets, double sigma_gen){
    AllJets.insert( AllJets.end(), jets.begin(), jets.end() );
    // di-jet asymmetry
    std::vector<Fjet> jjs;
    for (auto & J:jets){
        if (J.R<0.51 && J.R>0.49){
	    jjs.push_back(J);
	}
    }
    std::sort(jjs.begin(), jjs.end(), compare_jet);
    if (jjs.size()>=2){
        auto j1 = jjs[0], j2 = jjs[1];
	if (j1.pT>120. && j2.pT>50. && 
	    std::abs(j1.eta)<2. && std::abs(j2.eta)<2.){
	    double dphi = j1.phi - j2.phi;
	    if (std::cos(dphi) < std::cos(2./3.*M_PI)) {
                int index = 0;
		double x = std::abs(j1.pT-j2.pT)/(j1.pT+j2.pT);
		for (int i=0; i<xJbins.size()-1; i++){
	            if (xJbins[i]<x && x< xJbins[i+1]) {
	                index = i; 
			break;
		    }
		}
		xJ[index] += sigma_gen;
            }
	}
    }
    // shape
    for (auto & J:jets){
        if ((J.R<0.41) && (J.R>0.39) && (std::abs(J.eta) < 1.6) ) { // CMS cut
            int ii=0;
            for (int i=0; i<shape_NpT; i++){
                if ((shape_pTbins[i]<J.pT) && (J.pT<shape_pTbins[i+1])) {
                    ii = i;
                    break;
                }
            }
            if (ii==NpT) continue;
            for (int i=0; i<shape_Nr; i++){
                shapes[ii][i] += sigma_gen*J.shape[i];
                dnchdr[ii][i] += sigma_gen*J.dndr[i];
            }
            shape_w[ii] += sigma_gen;
        }
    }

    for (auto & J : jets){
    // check jet R
    int iR;
        for (int i=0; i<Rs.size(); i++){
            if ( ( (Rs[i]-.01)<J.R)
                && (J.R< (Rs[i]+.01))
               ) {
                iR = i;
                break;
            }
        }
        if (iR==Rs.size()) continue;

        if (std::abs(J.eta) < 2.1){
            // check flavor content
            bool isD = false;
            bool isDpTcut = false;
            bool isB = false;
            bool isBpTcut = false;
            for (auto & p : J.Ftags) {
                int f = std::abs(p.pid);
                isD = isD || ((f==411)||(f==421)||(f==413)||(f==423)||(f==4));
                isB = isB || ((f==511)||(f==521)||(f==513)||(f==523)||(f==5));
                isBpTcut = isBpTcut || (
                       ((f==511)||(f==521)||(f==513)||(f==523)||(f==5))
                            && (p.p.xT()>5.));
                isDpTcut = isDpTcut || (
                       ((f==411)||(f==421)||(f==413)||(f==423)||(f==4))
                         && (p.p.xT()>5.));
            }
            // check jet pT cut
            int ii=0;
            for (int i=0; i<NpT; i++){
                if ( (pTbins[i]<J.pT) && (J.pT<pTbins[i+1]) ) {
                    ii = i;
                    break;
                }
            }
            if (ii==NpT) continue;
            if (isB) dB0dpT[iR][ii] += sigma_gen;
            if (isBpTcut) dBdpT[iR][ii] += sigma_gen;
            if (isD && (!isB)) dD0dpT[iR][ii] += sigma_gen;
            if (isDpTcut && (!isBpTcut)) dDdpT[iR][ii] += sigma_gen;
	    if ((!isDpTcut)&&(!isBpTcut)) dsigmadpT[iR][ii] += sigma_gen;
        }
    }     
}

void JetStatistics::write(std::string fheader){
    for (int iR=0; iR<Rs.size(); iR++){
        std::stringstream filename;
        filename << fheader << "-jet-R-"<<iR<<"-"<<"spectra.dat";
        std::ofstream f(filename.str());
        for (int i=0; i<NpT; i++) {
            f    << (pTbins[i]+pTbins[i+1])/2. << " "
                 << pTbins[i] << " "
                 << pTbins[i+1] << " "
                 << dsigmadpT[iR][i]/binwidth[i] << " "
                 << dDdpT[iR][i]/binwidth[i] << " "
                 << dBdpT[iR][i]/binwidth[i] << " "
         << dD0dpT[iR][i]/binwidth[i] << " "
                 << dB0dpT[iR][i]/binwidth[i] << std::endl;
        }
        f.close();
    }

    std::stringstream filename2;
    filename2 << fheader << "-jetshape.dat";
    std::ofstream f2(filename2.str());
    f2 << "#";
    for (auto it : shape_rbins) f2 << it << " ";
    f2 << std::endl;
    for (int i=0; i<shape_NpT; i++) {
        f2   << (shape_pTbins[i]+shape_pTbins[i+1])/2. << " "
             << shape_pTbins[i] << " "
             << shape_pTbins[i+1] << " ";
    for (int j=0; j<shapes[i].size(); j++){
        f2 << shapes[i][j]/(shape_w[i]+1e-15) << " "; 
    }
    f2 << std::endl;
    }
    f2.close();

    std::stringstream filename3;
    filename3 << fheader << "-jetdnchdr.dat";
    std::ofstream f3(filename3.str());
    f3 << "#";
    for (auto it : shape_rbins) f3 << it << " ";
    f3 << std::endl;
    for (int i=0; i<shape_NpT; i++) {
        f3   << (shape_pTbins[i]+shape_pTbins[i+1])/2. << " "
             << shape_pTbins[i] << " "
             << shape_pTbins[i+1] << " ";
        for (int j=0; j<dnchdr[i].size(); j++){
            f3 << dnchdr[i][j]/(shape_w[i]+1e-15) << " ";
        }
        f3 << std::endl;
    }
    f3.close();

    std::stringstream filename4;
    filename4 << fheader << "-xJ.dat";
    std::ofstream f4(filename4.str());
    for (int i=0; i<xJbins.size()-1; i++) {
        f4   << (xJbins[i]+xJbins[i+1])/2. << " "
             << xJbins[i] << " "
             << xJbins[i+1] << " "
             << xJ[i] << std::endl;
    }
    f4.close();

    std::stringstream filename5;
    filename5 << fheader << "-allHF.dat";
    std::ofstream f5(filename5.str());
    for (auto & J: AllJets) {
	if (J.Ftags.size()==0) continue;
	double LpT = 0.; int index=-1;
	for (int i=0; i<J.Ftags.size(); i++) {
            auto p = J.Ftags[i];
	    if (p.p.xT()>LpT) {
		  LpT = p.p.xT();
		  index = i;
	    }
	}
	auto pf = J.Ftags[index];
        int f = std::abs(pf.pid), ftag;
        if ((f==411)||(f==421)||(f==413)||(f==423)||(f==4)) ftag=4;
        if ((f==511)||(f==521)||(f==513)||(f==523)||(f==5)) ftag=5;
        f5 << J.sigma << " " << J.R << " " 
           << J.pT << " " << J.eta << " " << J.phi << " " << J.M2  << " " 
           << f << " " << ftag << " " 
	   << pf.p.xT() << " " << pf.p.rap() << " " << pf.p.phi() << " ";
	for (auto & it : J.shape) {
	    f5 << it << " ";
	}
	f5 << std::endl;
    }
    f5.close();

}

