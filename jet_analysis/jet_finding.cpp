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
 MR("Gmu"),
 BRs({{411, 0.1655}, {421, 0.0666}, {431, 0.0867},
      {511, 0.2311}, {521, 0.2040}, {531, 0.1918}
     }){
    if (need_response_table) MR.init(table_path);
    else MR.load(table_path);

    HF2mu.readString("Random:setSeed = on");
    HF2mu.readString("Random:seed = 0");
    HF2mu.readString("ProcessLevel:all = off");
    HF2mu.readString("HadronLevel:all = on");
    HF2mu.readString("HadronLevel:decay = on");
    HF2mu.init();

}

JetFinder::JetFinder(int _Neta, int _Nphi, double _etamax)
:Neta(_Neta), Nphi(_Nphi), 
 etamax(_etamax), etamin(-_etamax), 
 phimax(M_PI), phimin(-M_PI), 
 deta((etamax-etamin)/(_Neta-1)), 
 dphi((phimax-phimin)/Nphi),
 MR("Gmu"),
 BRs({{411, 0.1655}, {421, 0.0666}, {431, 0.0867},
      {511, 0.2311}, {521, 0.2040}, {531, 0.1918}
     }){

    HF2mu.readString("Random:setSeed = on");
    HF2mu.readString("Random:seed = 0");
    HF2mu.readString("ProcessLevel:all = off");
    HF2mu.readString("HadronLevel:all = on");
    HF2mu.readString("HadronLevel:decay = on");
    HF2mu.init();
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
	    auto absid = std::abs(p.pid);
	    if (p.charged || (p.pid==421) ) {
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



void JetFinder::FindJets(std::vector<double> Rs_, 
                         double jetpTMin, 
                         double jetyMin, 
                         double jetyMax,
			 bool chg_trigger){
    Rs = Rs_;
    LOG_INFO << "find jet";
    /*LOG_INFO << "trigger by gamma";
    double phigamma = 100;
    for (auto & p: plist){
      if (p.T0<-10) {
	      phigamma = p.p.phi();
		      break;
      }
    }
    if (phigamma > 10) {
	    LOG_INFO << "WORONGGGGG!";
	    Jets.clear();
	    return;
    }*/
    Jets.clear();
    HFs.clear();
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
        // recombine all the bins with negative contribution, subtract the paddle
        for (auto & j : jets){
            Fjet J;
            J.sigma = sigma;
            fourvec jetP{j.e(), j.px(), j.py(), j.pz()};
            fourvec newP{0., 0., 0., 0.};
              double Jphi = jetP.phi();
	      //LOG_INFO << Jphi << " +++" << std::endl;
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
		// Further trigger on gamma-away side:
		//if (std::cos(J.phi-phigamma)<0.)
        	Jets.push_back(J);
	    }
	    else{
	        // trigger on high-pT constituents
	        /*bool trigger = false;
	        for (auto & p: plist) {
                    double absdphi = std::abs(Jphi - p.p.phi());
                    if (absdphi>M_PI) absdphi = 2.*M_PI-absdphi;
                    double deta = Jeta-p.p.pseudorap();
                    double dR = std::sqrt(absdphi*absdphi + deta*deta);
                    if ( (dR<J.R) && p.charged && p.p.xT()>5.) trigger = true;
	        }*/
	        //if (trigger) 
		Jets.push_back(J);   
	    }
        }
    }
}


void JetFinder::FindHF(std::vector<particle> plist) {
    // find heavy meson, decay it until it produce a muon (weighted by branching ratio)
    // the candidate contains the HF meson and its decay product
    for (auto & pIn : plist){
	int absid = std::abs(pIn.pid);
        bool isHF = (pIn.pid==421);
        if (isHF) HFs.push_back(pIn);
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
        // response effect to FF
	/*for(double eta=Jeta-J.R; eta<=Jeta+J.R; eta+=.1){
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
        }*/
        for (int i=0; i<J.dNdz.size(); i++)
            J.dNdz[i] /= (zbins[i+1]-zbins[i]);
        for (int i=0; i<J.dNdpT.size(); i++)
            J.dNdpT[i] /= (zpTbins[i+1]-zpTbins[i]);
    }

    // for HF FF
    // Remerber that jets are ordered from high to low pT
    // Only assign HF particles to higher-pT jet if
    // there is a conflict
    // Label all HFs as visuable
    for (int iR=0; iR<Rs.size(); iR++){
	double thisR = Rs[iR];
        for(auto & p : HFs) p.is_virtual = false;
        for (int i=0; i<Jets.size(); i++){
            auto & J = Jets[i];
       	    if ( (J.R < thisR-0.01) || (thisR+.01 < J.R) ) continue;
            double Jphi = J.phi;
            double Jeta = J.eta;
            J.dDdz.resize(zbins.size()-1);
            J.dBdz.resize(zbins.size()-1);
            for (auto & it : J.dDdz) it = 0.;
            for (auto & it : J.dBdz) it = 0.;
            for (auto & hf : HFs) {
                if (hf.is_virtual) continue;
            
                double absdphi = std::abs(Jphi - hf.p.phi());
                if (absdphi>M_PI) absdphi = 2.*M_PI-absdphi;
                double deta = Jeta-hf.p.pseudorap();
                double dR_hf_J = std::sqrt(absdphi*absdphi + deta*deta);

                int absid = std::abs(hf.pid);
                bool isD = (absid/100)%100 == 4;
                bool isB = (absid/100)%100 == 5;
                if ((dR_hf_J < thisR) && isD && J.flavor==4) {
                   
                
        	    hf.is_virtual = true;
                    //double z = (hf.p.x()*J.pmu.x()+hf.p.y()*J.pmu.y()+hf.p.z()*J.pmu.z())/(J.pmu.x()*J.pmu.x()+J.pmu.y()*J.pmu.y()+J.pmu.z()*J.pmu.z());
                    double z = hf.p.xT()*std::cos(dR_hf_J)/J.pT;
		    int index = -1;
                    for (int i=0; i<J.dDdz.size(); i++){
                        if (z>zbins[i] && z<=zbins[i+1]) {
                            index = i;
                            break;
                        }
                    }
                    
                    // some exp cuts of ALICE:
                    bool cuts = (5<J.pT && J.pT<7 && 2<hf.p.xT() && hf.p.xT() < 7) ||
                                (7<J.pT && J.pT<10 && 3<hf.p.xT() && hf.p.xT() < 10) || 
                                (10<J.pT && J.pT<15 && 5<hf.p.xT() && hf.p.xT() < 15) || 
                                (15<J.pT && J.pT<50 && 5<hf.p.xT() && hf.p.xT() < 36) ;
                    
                    
                    if (index >= 0 && cuts) J.dDdz[index] += 1.;
                }
                if ((dR_hf_J < thisR) && isB && J.flavor==5) {
	            hf.is_virtual = true;
                    double z = hf.p.xT()*std::cos(dR_hf_J)/J.pT;
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
}


void JetFinder::LabelFlavor(){
    // label all HF as real
    for (int i=0; i<Rs.size(); i++){
	    double thisR = Rs[i];
    for (auto & p: HFs) p.is_virtual = false;
    for (auto & J : Jets){
        if ( (J.R < thisR-0.01) || (thisR+.01 < J.R) ) continue;
	// jets are ordered from high to low pT
	// we only assign the heavy quark the high-pT
	// jet if there is a conflict
        double Jphi = J.phi;
        double Jeta = J.eta;
        bool isD = false;
        bool isB = false;
        for (auto & hf : HFs){
            double absdphi = std::abs(Jphi - hf.p.phi());
            if (absdphi>M_PI) absdphi = 2.*M_PI-absdphi;
            double deta = Jeta-hf.p.pseudorap();
            double dR_hf_J = std::sqrt(absdphi*absdphi + deta*deta); 

            if (dR_hf_J<thisR && (!hf.is_virtual) ) {
                int absid = std::abs(hf.pid);
                isD = isD || ( hf.pid==421 && 3.<hf.p.xT() && hf.p.xT()<36. ); 
                //isB = isB || (absid==511 || absid==521 || absid==513 || absid==523 || absid==533 || absid==531); 
		// if this HF already belongs to a higher-pT jet
		// do not make it visuable to lower-pT jet
		hf.is_virtual = true;
		break;
            }
        }
        //if (isB) J.flavor = 5;
        //else if 
        if (isD) J.flavor = 4;
        else J.flavor = -1;
    }
  }
}

LeadingParton::LeadingParton():
pTbins({1,2,3,4,5,6,8,10,12,16,20,24,30,40,50,60,80,100,120,150,200,300,500}){
    NpT = pTbins.size()-1;
    nchg.resize(NpT); for (auto & it:nchg) it=0.;
    npi.resize(NpT); for (auto & it:npi) it=0.;
    nstrange.resize(NpT); for (auto & it:nstrange) it=0.;
    nD.resize(NpT); for (auto & it:nD) it=0.;
    nB.resize(NpT); for (auto & it:nB) it=0.;
    v2chg.resize(NpT); 
    for (auto & it:v2chg) {
        it.resize(2);
        for (auto & iit : it) iit=0.;
    }
    v3chg.resize(NpT); 
    for (auto & it:v3chg) {
        it.resize(2);
        for (auto & iit : it) iit=0.;
    }
    v2pi.resize(NpT); 
    for (auto & it:v2pi) {
        it.resize(2);
        for (auto & iit : it) iit=0.;
    }
    v3pi.resize(NpT); 
    for (auto & it:v3pi) {
        it.resize(2);
        for (auto & iit : it) iit=0.;
    }
    v2strange.resize(NpT); 
    for (auto & it:v2strange) {
        it.resize(2);
        for (auto & iit : it) iit=0.;
    }
    v3strange.resize(NpT); 
    for (auto & it:v3strange) {
        it.resize(2);
        for (auto & iit : it) iit=0.;
    }
    v2D.resize(NpT); 
    for (auto & it:v2D) {
        it.resize(2);
        for (auto & iit : it) iit=0.;
    }
    v3D.resize(NpT); 
    for (auto & it:v3D) {
        it.resize(2);
        for (auto & iit : it) iit=0.;
    }
    v2B.resize(NpT); 
    for (auto & it:v2B) {
        it.resize(2);
        for (auto & iit : it) iit=0.;
    }
    v3B.resize(NpT); 
    for (auto & it:v3B) {
        it.resize(2);
        for (auto & iit : it) iit=0.;
    }

    binwidth.resize(NpT); 
    for (int i=0; i<NpT; i++) 
        binwidth[i]=pTbins[i+1]-pTbins[i];
}
void LeadingParton::add_event(std::vector<particle> plist, 
                            double sigma_gen){
    for (auto & p : plist){
        {
            int pid = std::abs(p.pid);
            if (pid<6 || pid==21){ 
                LOG_INFO << "Parton in FS: "<< pid << " " 
                         << p.p.xT() << " " << p.p.rap();
		continue;
	    }
            double pT = p.p.xT();

	    bool ispi = (pid == 211) && (std::abs(p.p.pseudorap()) < 1.) ;
	    bool isstrange = (pid == 321) && (std::abs(p.p.pseudorap()) < 1.) ;
            bool ischg = (pid==211 || pid==321 || pid==2212) && (std::abs(p.p.pseudorap()) < 1.) ;
            bool isD = (pid==411 || pid == 421 || pid == 413 || pid == 423) && (std::abs(p.p.rap()) < 1.) ;
            bool isB = (pid==511 || pid == 521 || pid == 513  || pid == 523) && (std::abs(p.p.rap()) < 2.4) ;
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
                v2pi[ii][0] += sigma_gen*std::cos(2*p.p.phi());
                v2pi[ii][1] += sigma_gen*std::sin(2*p.p.phi());
                v3pi[ii][0] += sigma_gen*std::cos(3*p.p.phi());
                v3pi[ii][1] += sigma_gen*std::sin(3*p.p.phi());
            }
            if (ischg){
                nchg[ii] += sigma_gen;
                v2chg[ii][0] += sigma_gen*std::cos(2*p.p.phi());
                v2chg[ii][1] += sigma_gen*std::sin(2*p.p.phi());
                v3chg[ii][0] += sigma_gen*std::cos(3*p.p.phi());
                v3chg[ii][1] += sigma_gen*std::sin(3*p.p.phi());
            }
            if (isstrange){
                nstrange[ii] += sigma_gen;
                v2strange[ii][0] += sigma_gen*std::cos(2*p.p.phi());
                v2strange[ii][1] += sigma_gen*std::sin(2*p.p.phi());
                v3strange[ii][0] += sigma_gen*std::cos(3*p.p.phi());
                v3strange[ii][1] += sigma_gen*std::sin(3*p.p.phi());
            }
            if (isD) {
                nD[ii] += sigma_gen;
                v2D[ii][0] += sigma_gen*std::cos(2*p.p.phi());
                v2D[ii][1] += sigma_gen*std::sin(2*p.p.phi());
                v3D[ii][0] += sigma_gen*std::cos(3*p.p.phi());
                v3D[ii][1] += sigma_gen*std::sin(3*p.p.phi());
            }
            if (isB) {
                nB[ii] += sigma_gen;
                v2B[ii][0] += sigma_gen*std::cos(2*p.p.phi());
                v2B[ii][1] += sigma_gen*std::sin(2*p.p.phi());
                v3B[ii][0] += sigma_gen*std::cos(3*p.p.phi());
                v3B[ii][1] += sigma_gen*std::sin(3*p.p.phi());
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
             << nchg[i]/binwidth[i] << " "
             << npi[i]/binwidth[i] << " "
             << nD[i]/binwidth[i] << " "
             << nB[i]/binwidth[i] << std::endl;
    }
    f.close();

    /*std::stringstream filename2;
    filename2 << fheader << "-LeadingQn.dat";
    std::ofstream f2(filename2.str());

    for (int i=0; i<NpT; i++) {
        f2   << (pTbins[i]+pTbins[i+1])/2. << " "
             << pTbins[i] << " "
             << pTbins[i+1] << " "
             << nchg[i] << " " 
             << v2chg[i][0] << " " << v2chg[i][1] << " "
             << v3chg[i][0] << " " << v3chg[i][1] << " "
             << npi[i] << " " 
             << v2pi[i][0] << " " << v2pi[i][1] << " "
             << v3pi[i][0] << " " << v3pi[i][1] << " "
             << nstrange[i] << " "
             << v2strange[i][0] << " " << v2strange[i][1] << " "
             << v3strange[i][0] << " " << v3strange[i][1] << " "
             << nD[i] << " " 
             << v2D[i][0] << " " << v2D[i][1] << " "
             << v3D[i][0] << " " << v3D[i][1] << " "
             << nB[i] << " " 
             << v2B[i][0] << " " << v2B[i][1] << " "
             << v3B[i][0] << " " << v3B[i][1] << std::endl;
    }
    f2.close();*/
}

JetStatistics::JetStatistics(
     std::vector<double> _Rs, 
     std::vector<double> _shape_rbins,
     std::vector<double> _Frag_zbins,
     std::vector<double> _Frag_zpTbins):
Rs(_Rs),
//pTbins({5,7,9,11,13,15,18,21,24,27,30,35,40,45,50}),
pTbins({4,6,8,10,12,16,20,25,30,35,40,50,60,80,100,120,140,160,180,200,220}),
xJ_pTbins({100,112,126,141,
	   158,178,200,224,
	   251,282,316,398,562}),
shape_pTbins({10,20,30,40,60,80,100,120,2500}),
//Frag_pTbins({10,20,30,40,60,70,100,126,158,200,251,316,398,600,800}),
Frag_pTbins({5,7,10,15,50}),
shape_rbins(_shape_rbins),
Frag_zbins(_Frag_zbins),
Frag_zpTbins(_Frag_zpTbins),
xJbins({0, .05,  .1, .15,  .2,
         .25,  .3, .35,  .4, .45,
          .5, .55,  .6, .65,  .7,
         .75,  .8, .85,  .9, .95, 1.0})
{ 

    xJ_W.resize(xJ_pTbins.size()-1);
    for (auto & it: xJ_W) it = 0.;
    xJ.resize(xJ_pTbins.size()-1);
    for (auto & it: xJ) {
        it.resize(xJbins.size()-1);
        for (auto & iit : it) iit = 0.;
    }
    NpT = pTbins.size()-1;
    leading_yield.resize(NpT);
    subleading_yield.resize(NpT);
    for (auto & it : leading_yield) it=0.;
    for (auto & it : subleading_yield) it=0.;
    shape_NpT = shape_pTbins.size()-1;
    shape_Nr = shape_rbins.size()-1;

    Frags.resize(Frag_pTbins.size()-1);
    Frags_pT.resize(Frag_pTbins.size()-1);
    Frags_W.resize(Frag_pTbins.size()-1);
    Frags_D.resize(Frag_pTbins.size()-1);
    Frags_D_pT.resize(Frag_pTbins.size()-1);
    Frags_D_W.resize(Frag_pTbins.size()-1);
    Frags_B.resize(Frag_pTbins.size()-1);
    Frags_B_pT.resize(Frag_pTbins.size()-1);
    Frags_B_W.resize(Frag_pTbins.size()-1);

    D_in_jet_W.resize(Rs.size());
    for (auto & it: D_in_jet_W) {
        it.resize(Frag_pTbins.size()-1);
        for (auto & iit : it) iit = 0.;
    }
    
    B_in_jet_W.resize(Rs.size());
    for (auto & it: B_in_jet_W) {
        it.resize(Frag_pTbins.size()-1);
	for (auto & iit : it) iit = 0.;
    }

    D_in_jet.resize(Rs.size());
    for (auto & it: D_in_jet){
	it.resize(Frag_pTbins.size()-1);
        for (auto & iit : it) {
            iit.resize(Frag_zbins.size()-1);
	    for (auto & iiit : iit) iiit=0.;
	}
    }

    B_in_jet.resize(Rs.size());
    for (auto & it: B_in_jet){
        it.resize(Frag_pTbins.size()-1);
        for (auto & iit : it) {
            iit.resize(Frag_zbins.size()-1);
            for (auto & iiit : iit) iiit=0.;
        }
    }



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
    for (auto & it : Frags_D_pT) {
        it.resize(Frag_zpTbins.size()-1);
        for (auto & itt : it) itt =0.;
    }

    for (auto & it : Frags_B_W) it=0.;
    for (auto & it : Frags_B) {
        it.resize(Frag_zbins.size()-1);
        for (auto & itt : it) itt =0.;
    }
    for (auto & it : Frags_B_pT) {
        it.resize(Frag_zpTbins.size()-1);
        for (auto & itt : it) itt =0.;
    }


    Shape_W.resize(shape_NpT); for (auto & it: Shape_W) it=0.;
    shapes.resize(shape_NpT); 
    for (auto & it:shapes) {
        it.resize(_shape_rbins.size()-1);
        for (auto & iit : it ) iit=0.;
    }
    Shape_D_W.resize(shape_NpT); for (auto & it: Shape_D_W) it=0.;
    Dshapes.resize(shape_NpT); 
    for (auto & it:Dshapes) {
        it.resize(_shape_rbins.size()-1);
        for (auto & iit : it ) iit=0.;
    }
    Shape_B_W.resize(shape_NpT); for (auto & it: Shape_B_W) it=0.;
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
    JQ2.resize(Rs.size());
    for (auto & it:JQ2) {
        it.resize(NpT);
        for (auto & iit : it ){
            iit.resize(2);
            iit[0] = 0.; iit[1] = 0.;
        }
    }
    JQ3.resize(Rs.size());
    for (auto & it:JQ3) {
        it.resize(NpT);
        for (auto & iit : it ){
            iit.resize(2);
            iit[0] = 0.; iit[1] = 0.;
        }
    }
    JQ4.resize(Rs.size());
    for (auto & it:JQ4) {
        it.resize(NpT);
        for (auto & iit : it ){
            iit.resize(2);
            iit[0] = 0.; iit[1] = 0.;
        }
    }

    binwidth.resize(NpT);
    for (int i=0; i<NpT; i++)
        binwidth[i]=pTbins[i+1]-pTbins[i];
}
//std::ofstream  fdijet("dijet-lowpT.dat");
void JetStatistics::add_event(std::vector<Fjet> jets, double sigma_gen, fourvec x0){
    // di-jet asymmetry
    /*std::vector<Fjet> jjs;
    for (auto & J:jets){
        if (J.R<0.41 && J.R>0.39){
            jjs.push_back(J);
        }
    }
    std::sort(jjs.begin(), jjs.end(), compare_jet);

    if (jjs.size()>=2){
        auto j1 = jjs[0], j2 = jjs[1];
        if (std::abs(j1.eta)<2.1 && std::abs(j2.eta)<2.1) {
            double dphi = j1.phi - j2.phi;
            if (std::cos(dphi) < std::cos(7./8.*M_PI))
                fdijet << j1.pT << " " << j2.pT << " " << sigma_gen << std::endl;
	    double x = j2.pT/j1.pT;
            if ((std::cos(dphi) < std::cos(7./8.*M_PI)) && (x>0.32)) {
	        // back-to-back dijets
		int pTindex = -1;
                for (int i=0; i<xJ_pTbins.size()-1; i++){
                    if (xJ_pTbins[i]<j1.pT && j1.pT< xJ_pTbins[i+1]) {
                        pTindex = i;
                        break;
                    }
                }
		if (pTindex>=0) {
                    int rindex = -1;
                    for (int i=0; i<xJbins.size()-1; i++){
                        if (xJbins[i]<x && x < xJbins[i+1]) {
                            rindex = i; 
                            break;
                        }
                    }
		    if (rindex>=0){
                        xJ[pTindex][rindex] += sigma_gen;
                        xJ_W[pTindex] += sigma_gen;
                    }
		}
                // leading/subleading jet yield
		for (int i=0; i<pTbins.size()-1; i++){
                     if (pTbins[i]<j1.pT && j1.pT< pTbins[i+1])
		        leading_yield[i] += sigma_gen;
	             if (pTbins[i]<j2.pT && j2.pT< pTbins[i+1])
                        subleading_yield[i] += sigma_gen;
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
	    if (J.flavor==4) Shape_D_W[ii] += sigma_gen;
	    if (J.flavor==5) Shape_B_W[ii] += sigma_gen;
	    Shape_W[ii] += sigma_gen;
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
            for (int j=0; j<Frags[ii].size(); j++){
                if (J.flavor==4) Frags_D[ii][j] += sigma_gen*J.dNdz[j];
		if (J.flavor==5) Frags_B[ii][j] += sigma_gen*J.dNdz[j];
		Frags[ii][j] += sigma_gen*J.dNdz[j];
	    }
            for (int j=0; j<Frags_pT[ii].size(); j++){
                if (J.flavor==4) Frags_D_pT[ii][j] += sigma_gen*J.dNdpT[j];
		if (J.flavor==5) Frags_B_pT[ii][j] += sigma_gen*J.dNdpT[j];
		Frags_pT[ii][j] += sigma_gen*J.dNdpT[j];
	    }
        }
    }
*/

    // Flavor in jet
    for (int iR=0; iR<Rs.size(); iR++){
	double R = Rs[iR]; 
        for (auto & J : jets){
	    bool interested = (J.R<R+.01) && (J.R>R-.01);
            if ( !interested ) continue;
            int ii=-1;
            for (int i=0; i<Frag_pTbins.size()-1; i++){
                if ((Frag_pTbins[i]<J.pT) && (J.pT<Frag_pTbins[i+1])) {
                    ii = i;
                    break;
                }
            }
            if (ii<0) continue;
	    if (J.flavor==4) D_in_jet_W[iR][ii] += J.sigma;
	    if (J.flavor==5) B_in_jet_W[iR][ii] += J.sigma;
            for (int j=0; j<D_in_jet[iR][ii].size(); j++){
                if (J.flavor==4) D_in_jet[iR][ii][j] += J.sigma*J.dDdz[j];
		if (J.flavor==5) B_in_jet[iR][ii][j] += J.sigma*J.dBdz[j];
	    }
        }
    }


    // Yield and Vn
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
        if (std::abs(J.eta) < 0.9-J.R){
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
            JQ2[iR][ii][0] += sigma_gen * std::cos(2.*J.phi);
            JQ2[iR][ii][1] += sigma_gen * std::sin(2.*J.phi);
            JQ3[iR][ii][0] += sigma_gen * std::cos(3.*J.phi);
            JQ3[iR][ii][1] += sigma_gen * std::sin(3.*J.phi);
	    JQ4[iR][ii][0] += sigma_gen * std::cos(4.*J.phi);
            JQ4[iR][ii][1] += sigma_gen * std::sin(4.*J.phi);
        }
    }     
}

void JetStatistics::write(std::string fheader){
    for (int iR=0; iR<Rs.size(); iR++){
        std::stringstream filename;
        
        filename << fheader << "-spectra.dat";
        std::cout << filename.str() <<"!!!!!!!!!!!!"<<std::endl;
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

        /*std::stringstream filename2;
        filename2 << fheader << "-jet-R-"<<iR<<"-"<<"qn.dat";
        std::ofstream f2(filename2.str());
        for (int i=0; i<NpT; i++) {
            f2   << (pTbins[i]+pTbins[i+1])/2. << " "
                 << pTbins[i] << " "
                 << pTbins[i+1] << " "
                 << dsigmadpT[iR][i] << " " 
                 << JQ2[iR][i][0] << " " << JQ2[iR][i][1] << " " 
                 << JQ3[iR][i][0] << " " << JQ3[iR][i][1] << " "
		 << JQ4[iR][i][0] << " " << JQ4[iR][i][1] 
                 <<  std::endl;
        }
        f2.close();*/

    }

    /*std::stringstream filename1;
    filename1 << fheader << "-jet-dNdz.dat";
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
    filename11 << fheader << "-D-jet-dNdz.dat";
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

    std::stringstream filename111;
    filename111 << fheader << "-D-jet-dNdpT.dat";
    std::ofstream f111(filename111.str());
    f111 << "#";
    for (auto it : Frag_zpTbins) f111 << it << " ";
    f111 << std::endl;
    for (int i=0; i<Frags_pT.size(); i++) {
        f111   << (Frag_pTbins[i]+Frag_pTbins[i+1])/2. << " "
             << Frag_pTbins[i] << " "
             << Frag_pTbins[i+1] << " " << Frags_D_W[i] << " ";
        for (int j=0; j<Frags_D_pT[i].size(); j++){
            f111 << Frags_D_pT[i][j] << " " ;
        }
        f111 << std::endl;
    }
    f111.close();

    std::stringstream filename12;
    filename12 << fheader << "-B-jet-dNdz.dat";
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

    std::stringstream filename121;
    filename121 << fheader << "-B-jet-dNdpT.dat";
    std::ofstream f121(filename121.str());
    f121 << "#";
    for (auto it : Frag_zpTbins) f121 << it << " ";
    f121 << std::endl;
    for (int i=0; i<Frags_pT.size(); i++) {
        f121   << (Frag_pTbins[i]+Frag_pTbins[i+1])/2. << " "
             << Frag_pTbins[i] << " "
             << Frag_pTbins[i+1] << " " << Frags_B_W[i] << " ";
        for (int j=0; j<Frags_B_pT[i].size(); j++){
            f121 << Frags_B_pT[i][j] << " " ;
        }
        f121 << std::endl;
    }
    f121.close();


    std::stringstream filename2;
    filename2 << fheader << "-jetshape.dat";
    std::ofstream f2(filename2.str());
    f2 << "#";
    for (auto it : shape_rbins) f2 << it << " ";
    f2 << std::endl;
    for (int i=0; i<shape_NpT; i++) {
        f2   << (shape_pTbins[i]+shape_pTbins[i+1])/2. << " "
             << shape_pTbins[i] << " "
             << shape_pTbins[i+1] << " " << Shape_W[i] << " ";
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
             << shape_pTbins[i+1] << " " << Shape_D_W[i] << " ";
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
             << shape_pTbins[i+1] << " " << Shape_B_W[i] << " ";
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
           << xJ_pTbins[i+1] << " " << xJ_W[i] << " ";
        for (int j=0; j<xJ[i].size(); j++) {
            f5 << xJ[i][j] << " ";
        }
	f5 << std::endl;
    }
    f5.close();

    std::stringstream filename55;
    filename55 << fheader << "-DijetYield.dat";
    std::ofstream f55(filename55.str());
    for (int i=0; i<NpT; i++) {
            f55    << (pTbins[i]+pTbins[i+1])/2. << " "
                 << pTbins[i] << " "
                 << pTbins[i+1] << " "
                 << leading_yield[i]/binwidth[i] << " "
                 << subleading_yield[i]/binwidth[i] << std::endl;
    }
    f55.close();
   
*/
    for (int iR=0; iR<Rs.size(); iR++){
    /*std::stringstream filename100;
    filename100 << fheader << "-D-in-jet-iR-"<<iR<<".dat";
    std::ofstream f100(filename100.str());
    f100 << "#";
    for (auto it : Frag_zbins) f100 << it << " ";
    f100 << std::endl;
    for (int i=0; i<D_in_jet[iR].size(); i++) {
        f100   << (Frag_pTbins[i]+Frag_pTbins[i+1])/2. << " "
             << Frag_pTbins[i] << " "
             << Frag_pTbins[i+1] << " " << D_in_jet_W[iR][i] << " ";
        for (int j=0; j<D_in_jet[iR][i].size(); j++){
            f100 << D_in_jet[iR][i][j] << " " ;
        }
        f100 << std::endl;
    }
    f100.close(); */

    std::stringstream filename200;
    filename200 << fheader  << "-HFFF.dat";
    std::ofstream f200(filename200.str());
    f200 << "#";
    for (auto it : Frag_zbins) f200 << it << " ";
    f200 << std::endl;
    for (int i=0; i<D_in_jet[iR].size(); i++) {
        f200   << (Frag_pTbins[i]+Frag_pTbins[i+1])/2. << " "
             << Frag_pTbins[i] << " "
             << Frag_pTbins[i+1] << " " << D_in_jet_W[iR][i] << " ";
        for (int j=0; j<D_in_jet[iR][i].size(); j++){
            f200 << D_in_jet[iR][i][j] << " " ;
        }
        f200 << std::endl;
    }
    f200.close();
   }
}

JetHFCorr::JetHFCorr(std::vector<double> _rbins):
pTHFbins({4,20,2000}), 
rbins(_rbins){
    dDdr_W.resize(pTHFbins.size());
    for (auto & it: dDdr_W) it=0.;
    dBdr_W.resize(pTHFbins.size());
    for (auto & it: dBdr_W) it=0.;

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
                    int absid = std::abs(p.pid);
                    if (absid == 4 || 
                        absid == 411 || absid == 421 ||
                        absid == 413 || absid == 423) {
                        dDdr[pTindex][rindex] += sigma_gen;
                        dDdr_W[pTindex] += sigma_gen;
                    }    
                    if (absid == 5 ||
                        absid == 511 || absid == 521 ||
                        absid == 513 || absid == 523) {
                        dBdr[pTindex][rindex] += sigma_gen;
                        dBdr_W[pTindex] += sigma_gen;
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
             << pTHFbins[i+1] << " " << dDdr_W[i] << " ";
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
             << pTHFbins[i+1] << " " << dBdr_W[i] << " ";
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

