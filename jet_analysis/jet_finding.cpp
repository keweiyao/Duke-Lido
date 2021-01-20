#include <fstream>
#include <string>
#include <math.h>
#include <vector>
#include "jet_finding.h"
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

void JetFinder::MakeETower(double _vradial, 
                           double _Tfreeze, 
                           double pTmin,
                           std::vector<particle> _plist, 
                           std::vector<current> Sources,
                           int coarse_level,
			   bool charged_jet
			   ){
    vradial = _vradial;
    Tfreeze = _Tfreeze;
    plist.clear();
    plist = _plist;
    // Initialzie empty towers
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
	if (charged_jet){
	    if (p.charged) {
                Pmu[ieta][iphi] = Pmu[ieta][iphi] + p.p;
	        if (p.p.xT()>pTmin) PT[ieta][iphi] += p.p.xT();
	    }
	}
	else {
            Pmu[ieta][iphi] = Pmu[ieta][iphi] + p.p;
            if (p.p.xT()>pTmin) PT[ieta][iphi] += p.p.xT();
	}
    }
    if (Sources.size()==0) return;
    // put soft energy-momentum deposition into the towers
    clist.clear();
    clist.resize(coarseNeta);
    for (int ieta=0; ieta<coarseNeta; ieta++){
        clist[ieta].etas = etamin+ieta*coarsedeta;
        clist[ieta].p = azero;
    }
    for (auto & s: Sources){
        int i = corp_index(s.etas, etamin, etamax, coarsedeta, coarseNeta);
        if (i>=0) clist[i].p = clist[i].p + s.p;
    }
    for (int ieta=0; ieta<coarseNeta; ieta++) {
        double eta = etamin+ieta*coarsedeta;
        for (int iphi=0; iphi<coarseNphi; iphi++) {
            double phi = phimin+iphi*coarsedphi;
            double dpT = 0., dpTcut = 0.;
            for (auto & s: clist){
                double etas = s.etas;
                if (charged_jet){
		    dpT += 2./3.*MR.get_dpT_dydphi(eta-etas, phi, 
				                   s.p, vradial, 0.);
                    dpTcut += 2./3.*MR.get_dpT_dydphi(eta-etas, phi, 
				               s.p, vradial, pTmin/Tfreeze);
		}
		else{
                    dpT += MR.get_dpT_dydphi(eta-etas, phi, 
                                                   s.p, vradial, 0.);
                    dpTcut += MR.get_dpT_dydphi(eta-etas, phi, 
                                               s.p, vradial, pTmin/Tfreeze);
		}

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
                         double jetyMax,
			 bool chg_trigger) {
    LOG_INFO << "find jet";
    Jets.clear();
    HFs.clear();
    HFaxis.clear();
    std::vector<fastjet::PseudoJet> fjInputs;
    for (int ieta=0; ieta<Pmu.size(); ieta++) {
        for (int iphi=0; iphi<Pmu[ieta].size(); iphi++) {
            auto iit = Pmu[ieta][iphi];
            // when we do the first jet finding,
            // only use towers with positive contribution
            if (iit.t()>0.) {
                fastjet::PseudoJet fp(iit.x(), iit.y(), 
                                          iit.z(), iit.t());
                fjInputs.push_back(fp);
            }
        }
    }

    fastjet::Selector select_rapidity = fastjet::SelectorRapRange(jetyMin, jetyMax);
    for (int iR=0; iR<Rs.size(); iR++){
	double jetRadius = Rs[iR];
	fastjet::JetDefinition jetDef(fastjet::genkt_algorithm, jetRadius, -1);
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
	    if (!chg_trigger){
        	Jets.push_back(J);
	    }
	    else{
	        // trigger on high-pT constituents
	        bool trigger = false;
	        for (auto & p: plist) {
                    double absdphi = std::abs(Jphi - p.p.phi());
                    if (absdphi>M_PI) absdphi = 2.*M_PI-absdphi;
                    double deta = Jeta-p.p.pseudorap();
                    double dR = std::sqrt(absdphi*absdphi + deta*deta);
                    if ( (dR<J.R) && p.charged && p.p.xT()>5.) trigger = true;
	        }
	        if (trigger) Jets.push_back(J);   
	    }
        }
    }
}


void JetFinder::FindHF(std::vector<particle> plist) {
    HFs.clear();
    // Find all heavy flavor particles
    for (auto & p : plist){
        int absid = std::abs(p.pid);
        int meson_specie = (std::abs(p.pid)/100)%10;
        int baryon_specie = (std::abs(p.pid)/1000)%10;
        if (absid == 4 || meson_specie==4 || baryon_specie==4||
            absid == 5 || meson_specie==5 || baryon_specie==5){
            HFs.push_back(p);
        }
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

void JetFinder::Frag(std::vector<double> zbins, std::vector<double> zpTbins){
    for (int i=0; i<Jets.size(); i++){
	auto & J = Jets[i];
        double Jphi = J.phi;
        double Jeta = J.eta;
        J.dNdz.resize(zbins.size()-1);
        for (auto & it : J.dNdz) it = 0.;
        J.dNdpT.resize(zpTbins.size()-1);
        for (auto & it : J.dNdpT) it = 0.;

	for (auto & p : plist) {
	    double absdphi = std::abs(Jphi - p.p.phi());
            if (absdphi>M_PI) absdphi = 2.*M_PI-absdphi;
            double deta = Jeta-p.p.pseudorap();
            double dR = std::sqrt(absdphi*absdphi + deta*deta); 
	    // D(z)
            if (dR<J.R && p.charged) {
                double z = std::min(p.p.xT()*std::cos(dR)/J.pT,.999);
                int index = -1;
                for (int i=0; i<J.dNdz.size(); i++){
                    if (z>zbins[i] && z<=zbins[i+1]) {
                        index = i;
                        break;
                    }
                }
		if (index >= 0) J.dNdz[index] += 1.;
            }
	    // DpT
            if (dR<J.R && p.charged) {
                double pT = p.p.xT();
                int index = -1;
                for (int i=0; i<J.dNdpT.size(); i++){
                    if (pT>zpTbins[i] && pT<=zpTbins[i+1]) {
                        index = i;
                        break;
                    }
                }
                if (index >= 0) J.dNdpT[index] += 1.;
            } 
        }
        // response effect
	for(double eta=Jeta-J.R; eta<=Jeta+J.R; eta+=.1){
            double Rp = std::sqrt(1e-9+std::pow(J.R,2)-std::pow(Jeta-eta,2)); 
	    for (double phi=Jphi-Rp; phi<=Jphi+Rp; phi+=.1){
		double dR = std::sqrt(std::pow(Jeta-eta,2)
				+ std::pow(Jphi-phi,2));
                double newphi = phi;
                if (phi<-M_PI) newphi+=2*M_PI;
                if (phi>M_PI) newphi-=2*M_PI;
                for (auto & s: clist){
                    double etas = s.etas;
		    // dNdz
		    for (int i=0; i<J.dNdz.size(); i++){
                        double pTl = zbins[i]*J.pT/std::cos(dR);
			double pTh = zbins[i+1]*J.pT/std::cos(dR);
			if(pTh/Tfreeze<30.){
                        double Yield_l = MR.get_dpT_dydphi(eta-etas, newphi, 
				    s.p,vradial, pTl/Tfreeze);
                        double Yield_h = MR.get_dpT_dydphi(eta-etas, newphi,
                                    s.p,vradial, pTh/Tfreeze);
		  	J.dNdz[i] += 2./3.*(Yield_l-Yield_h)
				       /((pTh+pTl)/2.)
				       *.1*.1;
			}
		    }
		    // dNdpT
                    for (int i=0; i<J.dNdpT.size(); i++){
                        double pTl = zpTbins[i];
                        double pTh = zpTbins[i+1];
                        if(pTh/Tfreeze<30.){
                        double Yield_l = MR.get_dpT_dydphi(eta-etas, newphi,
                                    s.p,vradial, pTl/Tfreeze);
                        double Yield_h = MR.get_dpT_dydphi(eta-etas, newphi,
                                    s.p,vradial, pTh/Tfreeze);
                        J.dNdpT[i] += 2./3.*(Yield_l-Yield_h)
                                       *.1*.1;
                        }
                    }
                }
	    }
        }
        for (int i=0; i<J.dNdz.size(); i++)
            J.dNdz[i] /= (zbins[i+1]-zbins[i]);
        for (int i=0; i<J.dNdpT.size(); i++)
            J.dNdpT[i] /= (zpTbins[i+1]-zpTbins[i]);
    }

    // for HF
    for (int i=0; i<Jets.size(); i++){
        auto & J = Jets[i];
        double Jphi = J.phi;
        double Jeta = J.eta;
        J.dDdz.resize(zbins.size()-1);
	J.dBdz.resize(zbins.size()-1);
        for (auto & it : J.dDdz) it = 0.;
        for (auto & it : J.dBdz) it = 0.;
        for (auto & p : HFs) {
            int absid = std::abs(p.pid);
            bool isD = 
                   (absid == 4 ||
                    absid == 411 || absid == 421 ||
                    absid == 413 || absid == 423) ;
            bool isB = 
                   (absid == 5 ||
                    absid == 511 || absid == 521 ||
                    absid == 513 || absid == 523) ;
            double absdphi = std::abs(Jphi - p.p.phi());
            if (absdphi>M_PI) absdphi = 2.*M_PI-absdphi;
            double deta = Jeta-p.p.pseudorap();
            double dR = std::sqrt(absdphi*absdphi + deta*deta);
            if ((dR < J.R) && isD) {
                double z = std::min(p.p.xT()*std::cos(dR)/J.pT,.999);
		int index = -1;
                for (int i=0; i<J.dDdz.size(); i++){
                    if (z>zbins[i] && z<=zbins[i+1]) {
                        index = i;
                        break;
                    }
                }
                if (index >= 0) J.dDdz[index] += 1.;
            }
            if ((dR < J.R) && isB) {
                double z = std::min(p.p.xT()*std::cos(dR)/J.pT,.999);
                int index = -1;
                for (int i=0; i<J.dBdz.size(); i++){ 
                    if (z>zbins[i] && z<=zbins[i+1]) {
                        index = i;
                        break;
                    }
                }
                if (index >= 0) J.dBdz[index] += 1.;
            }

        }
        for (int i=0; i<J.dDdz.size(); i++)
            J.dDdz[i] /= (zbins[i+1]-zbins[i]);
        for (int i=0; i<J.dBdz.size(); i++)
            J.dBdz[i] /= (zbins[i+1]-zbins[i]);
    }

}


void JetFinder::LabelFlavor(){
    for (auto & J : Jets){
        double Jphi = J.phi;
        double Jeta = J.eta;
        bool isD = false;
        bool isB = false;
        for (auto & p : HFs){
            double absdphi = std::abs(Jphi - p.p.phi());
            if (absdphi>M_PI) absdphi = 2.*M_PI-absdphi;
            double deta = Jeta-p.p.pseudorap();
            double dR = std::sqrt(absdphi*absdphi + deta*deta); 
            if (dR<J.R) {
                int absid = std::abs(p.pid);
                isD = isD || 
                   ((absid == 4 || 
                    absid == 411 || absid == 421 ||
                    absid == 413 || absid == 423) ); 
                isB = isB || (
                   (absid == 5 || 
                    absid == 511 || absid == 521 ||
                    absid == 513 || absid == 523) ); 
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
    nstrange.resize(NpT); for (auto & it:nstrange) it=0.;
    nD.resize(NpT); for (auto & it:nD) it=0.;
    nB.resize(NpT); for (auto & it:nB) it=0.;
    v2chg.resize(NpT); for (auto & it:v2chg) it=0.;
    v2pi.resize(NpT); for (auto & it:v2pi) it=0.;
    v2strange.resize(NpT); for (auto & it:v2strange) it=0.;
    v2D.resize(NpT); for (auto & it:v2D) it=0.;
    v2B.resize(NpT); for (auto & it:v2B) it=0.;
    binwidth.resize(NpT); 
    for (int i=0; i<NpT; i++) 
        binwidth[i]=pTbins[i+1]-pTbins[i];
}
void LeadingParton::add_event(std::vector<particle> plist, double sigma_gen, double maxPT){
    for (auto & p : plist){
        if (std::abs(p.p.pseudorap()) < 1.) {
            int pid = std::abs(p.pid);
            if (pid<6 || pid==21){ 
                LOG_INFO << "Parton in FS: "<< pid << " " 
                         << p.p.xT() << " " << p.p.rap();
		continue;
	    }
            double pT = p.p.xT();

	    bool ispi = pid==111 || pid == 211;
	    bool isstrange = pid==311 || pid == 321;
            bool ischg = p.charged;
            bool isD = pid==411 || pid == 421 || pid == 413 
                    || pid == 423 || pid == 4;
            bool isB = pid==511 || pid == 521 || pid == 513 
		    || pid == 523 || pid ==5;
            int ii=0;
            for (int i=0; i<NpT; i++){
                if ( (pTbins[i]<pT) && (pT<pTbins[i+1]) ) {
                    ii = i;
                    break;
                }
            }
            if (ii==NpT) ii=NpT-1;
            if (ispi) {
                npi[ii] += sigma_gen;
                v2pi[ii] += sigma_gen*std::cos(2*p.p.phi());
            }
            if (ischg){
                nchg[ii] += sigma_gen;
                v2chg[ii] += sigma_gen*std::cos(2*p.p.phi());
            }
            if (isstrange){
                nstrange[ii] += sigma_gen;
                v2strange[ii] += sigma_gen*std::cos(2*p.p.phi());
            }
            if (isD) {
                nD[ii] += sigma_gen;
                v2D[ii] += sigma_gen*std::cos(2*p.p.phi());
            }
            if (isB) {
                nB[ii] += sigma_gen;
                v2B[ii] += sigma_gen*std::cos(2*p.p.phi());
            }
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
             << nchg[i]/binwidth[i] << " " << v2chg[i] << " "
             << npi[i]/binwidth[i] << " " << v2pi[i] << " "
             << nstrange[i]/binwidth[i] << " "<< v2strange[i] << " "
             << nD[i]/binwidth[i] << " "<< v2D[i] << " "
             << nB[i]/binwidth[i] << " "<< v2B[i] << std::endl;
    }
    f.close();
}

JetStatistics::JetStatistics(
     std::vector<double> _pTbins, 
     std::vector<double> _Rs,
     std::vector<double> _shape_pTbins, 
     std::vector<double> _shape_rbins,
     std::vector<double> _Frag_pTbins,
     std::vector<double> _Frag_zbins,
     std::vector<double> _Frag_zpTbins,
     std::vector<double> _xJ_pTbins){ 
    DijetInfo.clear();
    xJ_pTbins = _xJ_pTbins;
    for (int i=0; i<21; i++) xJbins.push_back(i*0.05);
    xJ.resize(xJ_pTbins.size()-1);
    for (auto & it: xJ) {
        it.resize(xJbins.size()-1);
        for (auto & iit : it) iit = 0.;
    }
    pTbins = _pTbins;
    Rs = _Rs;
    shape_pTbins = _shape_pTbins;
    shape_rbins = _shape_rbins;
    NpT = pTbins.size()-1;
    shape_NpT = shape_pTbins.size()-1;
    shape_Nr = shape_rbins.size()-1;

    Frag_pTbins = _Frag_pTbins;
    Frag_zbins = _Frag_zbins;
    Frag_zpTbins = _Frag_zpTbins;

    Frags.resize(Frag_pTbins.size()-1);
    Frags_pT.resize(Frag_pTbins.size()-1);
    Frags_W.resize(Frag_pTbins.size()-1);
    Frags_D.resize(Frag_pTbins.size()-1);
    Frags_D_W.resize(Frag_pTbins.size()-1);
    Frags_B.resize(Frag_pTbins.size()-1);
    Frags_B_W.resize(Frag_pTbins.size()-1);

    for (auto & it : Frags_W) it=0.;
    for (auto & it : Frags) {
        it.resize(Frag_zbins.size()-1);
	for (auto & itt : it) itt =0.;
    }
    for (auto & it : Frags_pT) {
        it.resize(Frag_zpTbins.size()-1);
        for (auto & itt : it) itt =0.;
    }

    for (auto & it : Frags_D_W) it=0.;
    for (auto & it : Frags_D) {
        it.resize(Frag_zbins.size()-1);
        for (auto & itt : it) itt =0.;
    }

    for (auto & it : Frags_B_W) it=0.;
    for (auto & it : Frags_B) {
        it.resize(Frag_zbins.size()-1);
        for (auto & itt : it) itt =0.;
    }

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
void JetStatistics::add_event(std::vector<Fjet> jets, double sigma_gen, fourvec x0){
    // di-jet asymmetry
    std::vector<Fjet> jjs;
    for (auto & J:jets){
        if (J.R<0.41 && J.R>0.39){
            jjs.push_back(J);
        }
    }
    std::sort(jjs.begin(), jjs.end(), compare_jet);

    if (jjs.size()>=2){
        auto j1 = jjs[0], j2 = jjs[1];
        if (j1.pT>78 && j2.pT>34 &&
            std::abs(j1.eta)<2.1 && std::abs(j2.eta)<2.1) {
            double dphi = j1.phi - j2.phi;

            if (std::cos(dphi) < std::cos(7./8.*M_PI)) {
		    std::vector<double> entry{j1.pT, j1.phi, j1.eta,j2.pT,j2.phi, j2.eta,x0.x(),x0.y(),sigma_gen};
                DijetInfo.push_back(entry);
		int pTindex = -1;
                for (int i=0; i<xJ_pTbins.size()-1; i++){
                    if (xJ_pTbins[i]<j1.pT && j1.pT< xJ_pTbins[i+1]) {
                        pTindex = i;
                        break;
                    }
                }
		if (pTindex>=0) {
                    int rindex = -1;
                    double x = j2.pT/j1.pT;
                    for (int i=0; i<xJbins.size()-1; i++){
                        if (xJbins[i]<x && x< xJbins[i+1]) {
                            rindex = i; 
                            break;
                        }
                    }
		    if (rindex>=0)
                        xJ[pTindex][rindex] += sigma_gen;
		}
            }
        }
    }
    // shape
    for (auto & J:jets){
        if ((J.R<0.41) && (J.R>0.39) && (std::abs(J.eta) < 2.0) ) {
            int ii=-1;
            for (int i=0; i<shape_NpT; i++){
                if ((shape_pTbins[i]<J.pT) && (J.pT<shape_pTbins[i+1])) {
                    ii = i;
                    break;
                }
            }
            if (ii<0) continue;
            for (int i=0; i<shape_Nr; i++){
                shapes[ii][i] += sigma_gen*J.shape[i];
                if (J.flavor==4) Dshapes[ii][i] += sigma_gen*J.shape[i];
                if (J.flavor==5) Bshapes[ii][i] += sigma_gen*J.shape[i];
            }
        }
    }

    // Frag
    for (int i=0; i<jets.size(); i++){
	auto & J = jets[i];
	
        if ((J.R<0.41) && (J.R>0.39) && (std::abs(J.eta) < 2.1) ) {
            bool triggered = true;
            for (int j=0; j<jets.size(); j++){
                if (i==j) continue;
	        if (jets[j].R<0.39||jets[j].R>0.41) continue;
                auto & J2 = jets[j];
                double absdphi = std::abs(J.phi - J2.phi);
                if (absdphi>M_PI) absdphi = 2.*M_PI-absdphi;
                double deta = J.eta-J2.eta;
                double dR = std::sqrt(absdphi*absdphi + deta*deta);
                if ((dR<1.0) && (J2.pT>J.pT)) triggered = false;
            }
            if (!triggered) continue;	
            int ii=-1;
            for (int i=0; i<Frag_pTbins.size()-1; i++){
                if ((Frag_pTbins[i]<J.pT) && (J.pT<Frag_pTbins[i+1])) {
                    ii = i;
                    break;
                }
            }
            if (ii<0) continue;
	    if (J.flavor==4) Frags_D_W[ii] += sigma_gen;
	    if (J.flavor==5) Frags_B_W[ii] += sigma_gen;
	    Frags_W[ii] += sigma_gen;
            double sumb = 0.;
            double sumc = 0.;
            for (int j=0; j<Frags[ii].size(); j++){
                if (J.flavor==4) Frags_D[ii][j] += sigma_gen*J.dDdz[j];
		if (J.flavor==5) Frags_B[ii][j] += sigma_gen*J.dBdz[j];
		Frags[ii][j] += sigma_gen*J.dNdz[j];
		Frags_pT[ii][j] += sigma_gen*J.dNdpT[j];
	    }
        }
    }

    // Yield
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
            int ii=-1;
            for (int i=0; i<NpT; i++){
                if ( (pTbins[i]<J.pT) && (J.pT<pTbins[i+1]) ) {
                    ii = i;
                    break;
                }
            }
            if (ii<0) continue;
            if (J.flavor==4) dDdpT[iR][ii] += sigma_gen;
            if (J.flavor==5) dBdpT[iR][ii] += sigma_gen;
            dsigmadpT[iR][ii] += sigma_gen;
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

    std::stringstream filename1;
    filename1 << fheader << "-jetFrag.dat";
    std::ofstream f1(filename1.str());
    f1 << "#";
    for (auto it : Frag_zbins) f1 << it << " ";
    f1 << std::endl;
    for (int i=0; i<Frags.size(); i++) {
        f1   << (Frag_pTbins[i]+Frag_pTbins[i+1])/2. << " "
             << Frag_pTbins[i] << " "
             << Frag_pTbins[i+1] << " " << Frags_W[i] << " ";
        for (int j=0; j<Frags[i].size(); j++){
            f1 << Frags[i][j] << " " ;
        }
        f1 << std::endl;
    }
    f1.close();

    std::stringstream filename15;
    filename15 << fheader << "-jet-dNdpT.dat";
    std::ofstream f15(filename15.str());
    f15 << "#";
    for (auto it : Frag_zpTbins) f15 << it << " ";
    f15 << std::endl;
    for (int i=0; i<Frags_pT.size(); i++) {
        f15   << (Frag_pTbins[i]+Frag_pTbins[i+1])/2. << " "
             << Frag_pTbins[i] << " "
             << Frag_pTbins[i+1] << " " << Frags_W[i] << " ";
        for (int j=0; j<Frags_pT[i].size(); j++){
            f15 << Frags_pT[i][j] << " " ;
        }
        f15 << std::endl;
    }
    f15.close();


    std::stringstream filename11;
    filename11 << fheader << "-jetFragD.dat";
    std::ofstream f11(filename11.str());
    f11 << "#";
    for (auto it : Frag_zbins) f11 << it << " ";
    f11 << std::endl;
    for (int i=0; i<Frags_D.size(); i++) {
        f11   << (Frag_pTbins[i]+Frag_pTbins[i+1])/2. << " "
             << Frag_pTbins[i] << " "
             << Frag_pTbins[i+1] << " " << Frags_D_W[i] << " ";;
        for (int j=0; j<Frags_D[i].size(); j++){
            f11 << Frags_D[i][j] << " " ;
        }
        f11 << std::endl;
    }
    f11.close();

    std::stringstream filename12;
    filename12 << fheader << "-jetFragB.dat";
    std::ofstream f12(filename12.str());
    f12 << "#";
    for (auto it : Frag_zbins) f12 << it << " ";
    f12 << std::endl;
    for (int i=0; i<Frags_B.size(); i++) {
        f12   << (Frag_pTbins[i]+Frag_pTbins[i+1])/2. << " "
             << Frag_pTbins[i] << " "
             << Frag_pTbins[i+1] << " "  << Frags_B_W[i] << " ";
        for (int j=0; j<Frags_B[i].size(); j++){
            f12 << Frags_B[i][j] << " ";
        }
        f12 << std::endl;
    }
    f12.close();

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
    for (int i=0; i<xJ.size(); i++){
	f5 << (xJ_pTbins[i]+xJ_pTbins[i+1])/2. << " "
	   << xJ_pTbins[i] << " "
           << xJ_pTbins[i+1] << " ";
        for (int j=0; j<xJ[i].size(); j++) {
            f5 << xJ[i][j] << " ";
        }
	f5 << std::endl;
    }
    f5.close();
    /*std::stringstream filename6;
       filename6 << fheader << "-dijet-x0.dat";
    std::ofstream f6(filename6.str());
    for (auto & it: DijetInfo){
	    for (auto & iit: it){
		    f6 << iit << " ";
	    }f6<<std::endl;
    }*/

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
        if (std::abs(J.eta)<1.6 && std::abs(J.pT)>60. 
	    && J.R<0.31 && J.R>0.29){
            for (auto & p : HFs) {
                if (std::abs(p.p.rap())<2.){
                    double absdphi = std::abs(J.phi - p.p.phi());
                    if (absdphi>M_PI) absdphi = 2.*M_PI-absdphi;
                    double deta = J.eta-p.p.pseudorap();
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
                    int absid = p.pid;
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

void TestSource(
             MediumResponse MR,
             std::vector<current> SourceList,
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
                double etas = s.etas;
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

