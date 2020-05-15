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


JetFinder::JetFinder(int _Neta, int _Nphi, double _etamax, bool need_response_table, std::string table_path)
:Neta(_Neta), Nphi(_Nphi), 
 etamax(_etamax), etamin(-_etamax), 
 phimax(M_PI), phimin(-M_PI), 
 deta((etamax-etamin)/(_Neta-1)), 
 dphi((phimax-phimin)/Nphi),
 MR("Gmu") {
    if (need_response_table) MR.init(table_path);
    else MR.load(table_path);    
}

JetFinder::JetFinder(int _Neta, int _Nphi, double _etamax)
:Neta(_Neta), Nphi(_Nphi), 
 etamax(_etamax), etamin(-_etamax), 
 phimax(M_PI), phimin(-M_PI), 
 deta((etamax-etamin)/(_Neta-1)), 
 dphi((phimax-phimin)/Nphi),
 MR("Gmu") {
}

void JetFinder::MakeETower(double vradial, 
                           double Tfreeze, 
                           double pTmin,
                           std::vector<particle> plist, 
                           std::vector<current> TypeOneSources,
                           std::vector<HadronizeCurrent> TypeTwoSources,
                           int coarse_level){
    // Initialzie empty towers
    LOG_INFO << "do jet finding";
    PT.clear();
    PT.resize(Neta); 
    for (auto & it : PT) {
        it.resize(Nphi);
        for (auto & iit: it) iit = 0.;
    }    
    Pmu.clear();
    Pmu.resize(Neta); 
    fourvec azero{0., 0., 0., 0.};
    for (auto & it : Pmu) {
        it.resize(Nphi);
        for (auto & iit: it) iit = azero;
    }    
    // Initialize coarse towerse for medium response contribution
    int coarseNeta = int(Neta/coarse_level), 
        coarseNphi = int(Nphi/coarse_level);    
    double coarsedeta = (etamax-etamin)/(coarseNeta-1); 
    double coarsedphi = (phimax-phimin)/coarseNphi;
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
        Pmu[ieta][iphi] = Pmu[ieta][iphi] + p.p;
        if (p.p.xT()>pTmin) PT[ieta][iphi] += p.p.xT();
    }
    // put soft energy-momentum deposition into the towers
    std::vector<current> clist;
    clist.resize(coarseNeta);
    for (int ieta=0; ieta<coarseNeta; ieta++){
        clist[ieta].chetas = std::cosh(etamin+ieta*coarsedeta);
        clist[ieta].shetas = std::sinh(etamin+ieta*coarsedeta);
        clist[ieta].p = azero;
    }
    for (auto & s: TypeOneSources){
        int i = corp_index(std::atanh(s.shetas/s.chetas), etamin, etamax, coarsedeta, coarseNeta);
        if (i>=0) clist[i].p = clist[i].p + s.p;
    }
    for (int ieta=0; ieta<coarseNeta; ieta++) {
        double eta = etamin+ieta*coarsedeta;
        for (int iphi=0; iphi<coarseNphi; iphi++) {
            double phi = phimin+iphi*coarsedphi;
            double dpT = 0., dpTcut = 0.;
            for (auto & s: clist){
                double etas = std::atanh(s.shetas/s.chetas);
                dpT += MR.get_dpT_dydphi(eta-etas, phi, s.p, vradial, 0.);
                dpTcut += MR.get_dpT_dydphi(eta-etas, phi, s.p, vradial, pTmin/Tfreeze);
            }
            coarsePT[0][ieta][iphi] = dpT;
            coarsePT[1][ieta][iphi] = dpTcut;
        }
    }
    // Interpolate the coarse grid into the finer grid for jet finding
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
                    Pmu[ieta][iphi] = Pmu[ieta][iphi] + 
                                      Nmu0*coarsePT[0][ii+k1][jj+k2]*deta*dphi*u[k1]*v[k2];
                    PT[ieta][iphi] += coarsePT[1][ii+k1][jj+k2]*deta*dphi*u[k1]*v[k2];
                }
            }
        }
    }
    // done
}



void JetFinder::FindJets(std::vector<double> Rs, 
                         double jetpTMin, 
                         double jetyMin, 
                         double jetyMax) {
    Jets.clear();
    HFs.clear();
    HFaxis.clear();
    // Use the towers to do jet finding: anti-kT
    int power = -1; 
    // loop over jetradius
    for (int iR=0; iR<Rs.size(); iR++){
        double jetRadius = Rs[iR];
        fastjet::JetDefinition jetDef(fastjet::genkt_algorithm, jetRadius, power);
        std::vector<fastjet::PseudoJet> fjInputs;
        fastjet::Selector select_rapidity = fastjet::SelectorRapRange(jetyMin, jetyMax);
        for (auto & it : Pmu) {
            for (auto & iit : it) {
                // when we do the first jet finding,
                // only use towers with positive contribution
                if (iit.t()>0.){
                    fastjet::PseudoJet fp(iit.x(), iit.y(), 
                                          iit.z(), iit.t());
                    fp.set_user_info(new MyInfo(0, 0, 0));
                    fjInputs.push_back(fp);
                }
            }
        }
        fastjet::ClusterSequence clustSeq(fjInputs, jetDef);
        std::vector<fastjet::PseudoJet> AllJets = clustSeq.inclusive_jets(jetpTMin);
        std::vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(select_rapidity(AllJets));

        // once the jet is find, 
        // recombine all the bins with negative contribution
        for (auto & j : jets){
            Fjet J;
            J.sigma = sigma;
            fourvec jetP{j.e(), j.px(), j.py(), j.pz()};
            fourvec newP{0., 0., 0., 0.};
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
                    newP = newP + Pmu[ieta][iphi];
                }  
            } 
            J.pmu = newP;
            J.R = jetRadius;
            J.pT = newP.xT();
            J.phi = newP.phi();
            J.eta = newP.pseudorap();
            J.M2 = newP.m2();
            Jets.push_back(J);   
        }
    }
}


void JetFinder::FindHF(std::vector<particle> plist) {
    // Find all heavy flavor particles
    for (auto & p : plist){
        int absid = std::abs(p.pid);
        if (absid == 4 || 
            absid == 411 || absid == 421 ||
            absid == 413 || absid == 423 ||
            absid == 5 || 
            absid == 511 || absid == 521 ||
            absid == 513 || absid == 523) {
            HFs.push_back(p);
        }
    }
}
void JetFinder::CorrHFET(std::vector<double> rbins){
    // Correlate the energy production around each heavy flavor particles
    for (auto & p : HFs){
        int absid = std::abs(p.pid);
        Fjet J;
        if (absid == 4 || 
            absid == 411 || absid == 421 ||
            absid == 413 || absid == 423) J.flavor = 4;
        if (absid == 5 || 
            absid == 511 || absid == 521 ||
            absid == 513 || absid == 523) J.flavor = 5;
        J.R = 0.;
        J.pmu = p.p;
        J.pT = p.p.xT();
        J.eta = p.p.pseudorap();
        J.phi = p.p.phi();
        J.M2 = dot(p.p, p.p);
        J.sigma = sigma;
        fourvec p00{0,0,0,0};
        J.shape.resize(rbins.size()-1);
        for (auto & it : J.shape) it = 0.;
        for (double eta=J.eta-3.; eta<=J.eta+3.; eta+=deta){
            int ieta = corp_index(eta, etamin, etamax, deta, Neta);
            if (ieta<0) continue;
            double Rp = std::sqrt(std::pow(3.,2)-std::pow(J.eta-eta,2) ); 
            for (double phi=J.phi-Rp; phi<=J.phi+Rp; phi+=dphi){
                double dist = std::sqrt(std::pow(J.eta-eta,2)
                                      + std::pow(J.phi-phi,2));
                int index = -1;
                for (int i=0; i<J.shape.size(); i++){
                    if (dist>=rbins[i] && dist<rbins[i+1]) {
                        index = i;
                        break;
                    }
                }
                if (index==-1) continue;
                double newphi = phi;
                if (phi<-M_PI) newphi+=2*M_PI;
                if (phi>M_PI) newphi-=2*M_PI;
                int iphi = corp_index(newphi, phimin, phimax, dphi, Nphi);
                J.shape[index] += PT[ieta][iphi];
                if (dist<.1)
                     p00 = p00 + Pmu[ieta][iphi];
            }  
        }
        LOG_INFO << p.p.xT() << " vs " << p00.xT();
        LOG_INFO << p.p.pseudorap() << " " << p.p.phi() << " // "
                 << p00.pseudorap() << " " << p00.phi() ;
        // take out auto-correlation
        J.shape[0] -= J.pT; 
        for (int i=0; i<J.shape.size(); i++) 
            J.shape[i] /= (rbins[i+1]-rbins[i]);
        HFaxis.push_back(J);
    }
}

void JetFinder::CalcJetshape(std::vector<double> rbins){
    for (auto & J : Jets){
        double Jphi = J.phi;
        double Jeta = J.eta;
        J.shape.resize(rbins.size()-1);
        for (auto & it : J.shape) it = 0.;
        for (double eta=Jeta-3.; eta<=Jeta+3.; eta+=deta){
            int ieta = corp_index(eta, etamin, etamax, deta, Neta);
            if (ieta<0) continue;
            double Rp = std::sqrt(1e-9+std::pow(3.,2)-std::pow(Jeta-eta,2) ); 
            for (double phi=Jphi-Rp; phi<=Jphi+Rp; phi+=dphi){
                double dist = std::sqrt(std::pow(Jeta-eta,2) + std::pow(Jphi-phi,2));
                int index = -1;
                for (int i=0; i<J.shape.size(); i++){
                    if (dist>=rbins[i] && dist<rbins[i+1]) {
                        index = i;
                        break;
                    }
                }
                if (index==-1) continue;
                double newphi = phi;
                if (phi<-M_PI) newphi+=2*M_PI;
                if (phi>M_PI) newphi-=2*M_PI;
                int iphi = corp_index(newphi, phimin, phimax, dphi, Nphi);
                J.shape[index] += PT[ieta][iphi];
            }  
        }
        for (int i=0; i<J.shape.size(); i++) 
            J.shape[i] /= (rbins[i+1]-rbins[i]);
    }
}

void JetFinder::LabelFlavor(){
    for (auto & J : Jets){
        double Jphi = J.phi;
        double Jeta = J.eta;
        bool isD = false;
        bool isB = false;
        for (auto & p : HFs){
            double fphi = p.p.phi();
            double feta = p.p.rap();
            double absdphi = std::abs(Jphi - fphi);
            if (absdphi>M_PI) absdphi = 2.*M_PI-absdphi;
            double deta = Jeta-feta;
            double dR = std::sqrt(absdphi*absdphi + deta*deta);
            
            if (dR<0.3) {
                int absid = std::abs(p.pid);
                isD = isD || 
                   (absid == 4 || 
                    absid == 411 || absid == 421 ||
                    absid == 413 || absid == 423) 
                    && (p.p.xT() > 5.);
                isB = isB || 
                   (absid == 5 || 
                    absid == 511 || absid == 521 ||
                    absid == 513 || absid == 523) 
                    && (p.p.xT() > 5.);
            }
        }
        if (isB) J.flavor = 5;
        if (isD && (!isB)) J.flavor = 4;
        if ((!isD) && (!isB)) J.flavor = -1;
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
        if (std::abs(p.p.rap()) < 1.) {
            int pid = std::abs(p.pid);
            if (pid<6 || pid==21) 
                LOG_INFO << pid << " " 
                         << p.p.xT() << " " << p.p.rap();
            bool ispi = pid==111 || pid == 211;
            bool ischg = pid==111 || pid == 211 || pid == 311 
                 || pid== 321 || pid==2212;
            bool isD = pid==411 || pid == 421 || pid == 413 
                || pid == 423 || pid ==4;
            bool isB = pid==511 || pid == 521 || pid == 513 || pid == 523 || pid ==5;
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
JetStatistics::JetStatistics(
     std::vector<double> _pTbins, 
     std::vector<double> _Rs,
     std::vector<double> _shape_pTbins, 
     std::vector<double> _shape_rbins){
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

    shapes.resize(shape_NpT); 
    for (auto & it:shapes) {
        it.resize(_shape_rbins.size()-1);
        for (auto & iit : it ) iit=0.;
    }
    Dshapes.resize(shape_NpT); 
    for (auto & it:Dshapes) {
        it.resize(_shape_rbins.size()-1);
        for (auto & iit : it ) iit=0.;
    }
    Bshapes.resize(shape_NpT); 
    for (auto & it:Bshapes) {
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

    binwidth.resize(NpT);
    for (int i=0; i<NpT; i++)
        binwidth[i]=pTbins[i+1]-pTbins[i];
}
void JetStatistics::add_event(std::vector<Fjet> jets, double sigma_gen){
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
                if (J.flavor==4) Dshapes[ii][i] += sigma_gen*J.shape[i];
                if (J.flavor==5) Bshapes[ii][i] += sigma_gen*J.shape[i];
            }
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
            // check jet pT cut
            int ii=0;
            for (int i=0; i<NpT; i++){
                if ( (pTbins[i]<J.pT) && (J.pT<pTbins[i+1]) ) {
                    ii = i;
                    break;
                }
            }
            if (ii==NpT) continue;
            if (J.flavor==4) dDdpT[iR][ii] += sigma_gen;
            if (J.flavor==5) dBdpT[iR][ii] += sigma_gen;
            if (J.flavor==-1) dsigmadpT[iR][ii] += sigma_gen;
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
                 << dBdpT[iR][i]/binwidth[i] <<  std::endl;
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
            f2 << shapes[i][j] << " "; 
        }
        f2 << std::endl;
    }
    f2.close();

    std::stringstream filename3;
    filename3 << fheader << "-D-jetshape.dat";
    std::ofstream f3(filename3.str());
    f3 << "#";
    for (auto it : shape_rbins) f3 << it << " ";
    f3 << std::endl;
    for (int i=0; i<shape_NpT; i++) {
        f3   << (shape_pTbins[i]+shape_pTbins[i+1])/2. << " "
             << shape_pTbins[i] << " "
             << shape_pTbins[i+1] << " ";
        for (int j=0; j<Dshapes[i].size(); j++){
            f3 << Dshapes[i][j] << " "; 
        }
        f3 << std::endl;
    }
    f3.close();

    std::stringstream filename4;
    filename4 << fheader << "-B-jetshape.dat";
    std::ofstream f4(filename4.str());
    f4 << "#";
    for (auto it : shape_rbins) f4 << it << " ";
    f4 << std::endl;
    for (int i=0; i<shape_NpT; i++) {
        f4   << (shape_pTbins[i]+shape_pTbins[i+1])/2. << " "
             << shape_pTbins[i] << " "
             << shape_pTbins[i+1] << " ";
        for (int j=0; j<Bshapes[i].size(); j++){
            f4 << Bshapes[i][j] << " "; 
        }
        f4 << std::endl;
    }
    f4.close();

    std::stringstream filename5;
    filename5 << fheader << "-xJ.dat";
    std::ofstream f5(filename5.str());
    for (int i=0; i<xJbins.size()-1; i++) {
        f5   << (xJbins[i]+xJbins[i+1])/2. << " "
             << xJbins[i] << " "
             << xJbins[i+1] << " "
             << xJ[i] << std::endl;
    }
    f5.close();
}

JetHFCorr::JetHFCorr(std::vector<double> _pTHFbins, 
                     std::vector<double> _rbins)
:pTHFbins(_pTHFbins), rbins(_rbins){
    dDdr.clear();
    dDdr.resize(pTHFbins.size());
    for (auto & it: dDdr){
        it.resize(rbins.size());
        for (auto & iit: it) iit = 0.;
    }
    dBdr.clear();
    dBdr.resize(pTHFbins.size());
    for (auto & it: dBdr){
        it.resize(rbins.size());
        for (auto & iit: it) iit = 0.;
    }
}

void JetHFCorr::add_event(std::vector<Fjet> jets, 
                  std::vector<particle> HFs, 
                  double sigma_gen){
    for (auto & J : jets){
        if (std::abs(J.eta)<1.6 && std::abs(J.pT)>60. && J.R<0.31 && J.R>0.29){
            for (auto & p : HFs) {
                if (std::abs(p.p.rap())<2.){
                    double absdphi = std::abs(J.phi - p.p.phi());
                    if (absdphi>M_PI) absdphi = 2.*M_PI-absdphi;
                    double deta = J.eta-p.p.rap();
                    double dR = std::sqrt(absdphi*absdphi + deta*deta);
                    // find r index
                    int rindex = -1;
                    for (int i=0; i<rbins.size()-1; i++){
                        if (rbins[i]<dR && dR<rbins[i+1]){
                            rindex = i;
                            break;
                        }
                    }
                    if (rindex==-1) continue;
                    // find pT index
                    int pTindex = -1;
                    for (int i=0; i<pTHFbins.size()-1; i++){
                        if (pTHFbins[i]<p.p.xT() && p.p.xT()<pTHFbins[i+1]){
                            pTindex = i;
                            break;
                        }
                    }
                    if (pTindex==-1) continue;      
                    int absid = std::abs(p.pid);
                    if (absid == 4 || 
                        absid == 411 || absid == 421 ||
                        absid == 413 || absid == 423) {
                        dDdr[pTindex][rindex] += sigma_gen;
                    }    
                    if (absid == 5 || 
                        absid == 511 || absid == 521 ||
                        absid == 513 || absid == 523) {
                        dBdr[pTindex][rindex] += sigma_gen;
                    }          
                }
            }
        }
    }
}

void JetHFCorr::write(std::string fheader){
    std::stringstream filename1;
    filename1 << fheader << "-dDdr.dat";
    std::ofstream f1(filename1.str());
    f1 << "#";
    for (auto it : rbins) f1 << it << " ";
    f1 << std::endl;
    for (int i=0; i<pTHFbins.size()-1; i++) {
        f1   << (pTHFbins[i]+pTHFbins[i+1])/2. << " "
             << pTHFbins[i] << " "
             << pTHFbins[i+1] << " ";
        for (int j=0; j<rbins.size()-1; j++){
            f1 << dDdr[i][j] << " "; 
        }
        f1 << std::endl;
    }
    f1.close();

    std::stringstream filename2;
    filename2 << fheader << "-dBdr.dat";
    std::ofstream f2(filename2.str());
    f2 << "#";
    for (auto it : rbins) f2 << it << " ";
    f2 << std::endl;
    for (int i=0; i<pTHFbins.size()-1; i++) {
        f2   << (pTHFbins[i]+pTHFbins[i+1])/2. << " "
             << pTHFbins[i] << " "
             << pTHFbins[i+1] << " ";
        for (int j=0; j<rbins.size()-1; j++){
            f2 << dBdr[i][j] << " "; 
        }
        f2 << std::endl;
    }
    f2.close();
}


HFETCorr::HFETCorr(std::vector<double> _pTHFbins, 
                     std::vector<double> _rbins)
:pTHFbins(_pTHFbins), rbins(_rbins){
    D_dPTdr.clear();
    D_dPTdr.resize(pTHFbins.size());
    for (auto & it: D_dPTdr){
        it.resize(rbins.size());
        for (auto & iit: it) iit = 0.;
    }
    B_dPTdr.clear();
    B_dPTdr.resize(pTHFbins.size());
    for (auto & it: B_dPTdr){
        it.resize(rbins.size());
        for (auto & iit: it) iit = 0.;
    }
}

void HFETCorr::add_event(
                  std::vector<Fjet> HFaxis, 
                  double sigma_gen){
    for (auto & J : HFaxis){
        if (std::abs(J.eta)>2.1) continue;
        // find pT index
        int pTindex = -1;
        for (int i=0; i<pTHFbins.size()-1; i++){
            if (pTHFbins[i]<J.pT && J.pT<pTHFbins[i+1]){
                pTindex = i;
                break;
            }
        }
        if (pTindex==-1) continue;   
        if (J.flavor==4) {
            for (int j=0; j<rbins.size()-1; j++)
                D_dPTdr[pTindex][j] += sigma_gen*J.shape[j];
        }    
        if (J.flavor==5) {
            for (int j=0; j<rbins.size()-1; j++)
                B_dPTdr[pTindex][j] += sigma_gen*J.shape[j];
        } 
    }
}

void HFETCorr::write(std::string fheader){
    std::stringstream filename1;
    filename1 << fheader << "-D-dpTdr.dat";
    std::ofstream f1(filename1.str());
    f1 << "#";
    for (auto it : rbins) f1 << it << " ";
    f1 << std::endl;
    for (int i=0; i<pTHFbins.size()-1; i++) {
        f1   << (pTHFbins[i]+pTHFbins[i+1])/2. << " "
             << pTHFbins[i] << " "
             << pTHFbins[i+1] << " ";
        for (int j=0; j<rbins.size()-1; j++){
            f1 << D_dPTdr[i][j] << " "; 
        }
        f1 << std::endl;
    }
    f1.close();

    std::stringstream filename2;
    filename2 << fheader << "-B-dpTdr.dat";
    std::ofstream f2(filename2.str());
    f2 << "#";
    for (auto it : rbins) f2 << it << " ";
    f2 << std::endl;
    for (int i=0; i<pTHFbins.size()-1; i++) {
        f2   << (pTHFbins[i]+pTHFbins[i+1])/2. << " "
             << pTHFbins[i] << " "
             << pTHFbins[i+1] << " ";
        for (int j=0; j<rbins.size()-1; j++){
            f2 << B_dPTdr[i][j] << " "; 
        }
        f2 << std::endl;
    }
    f2.close();
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

