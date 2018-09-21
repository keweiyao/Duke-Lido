#include "workflow.h"
#include "simpleLogger.h"
#include <fstream>
#include "random.h"
#include "matrix_elements.h"
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/variant/get.hpp>
#include <boost/any.hpp>
#include <boost/foreach.hpp>
#include "logo.h"
#include "Langevin.h"

std::map<int, std::vector<Process> > AllProcesses;

void init_process(Process& r, std::string mode){
     switch(r.which()){
			case 0:
				if (boost::get<Rate22>(r).IsActive())
					if(mode == "new"){
						boost::get<Rate22>(r).initX("table.h5");
						boost::get<Rate22>(r).init("table.h5");
					} else{
						boost::get<Rate22>(r).loadX("table.h5");
						boost::get<Rate22>(r).load("table.h5");
					}
				else return;
				break;
			case 1:
				if (boost::get<Rate23>(r).IsActive())
					if(mode == "new"){
						boost::get<Rate23>(r).initX("table.h5");
						boost::get<Rate23>(r).init("table.h5");
					} else{
						boost::get<Rate23>(r).loadX("table.h5");
						boost::get<Rate23>(r).load("table.h5");
					}
				else return;
				break;
			case 2:
				if (boost::get<Rate32>(r).IsActive())
						if(mode == "new"){
								boost::get<Rate32>(r).initX("table.h5");
								boost::get<Rate32>(r).init("table.h5");
						} else{
								boost::get<Rate32>(r).loadX("table.h5");
								boost::get<Rate32>(r).load("table.h5");
						}
				else return;
				break;
			 case 3:
				if (boost::get<Rate12>(r).IsActive())
					if (mode == "new") {
						boost::get<Rate12>(r).init("table.h5");
				    } else{
						boost::get<Rate12>(r).load("table.h5");
				    }
				else return;
				break;
			 case 4:
				if (boost::get<Rate21>(r).IsActive())
					if (mode == "new") {
						boost::get<Rate21>(r).init("table.h5");
					} else {
						boost::get<Rate21>(r).load("table.h5");
					}
				else return;
				break;
			default:
				exit(-1);
				break;
		}
}


void initialize(std::string mode, std::string path, 
				double mu, double const_alpha_s, double A, double B){
	print_logo();
    initialize_mD_and_scale(0, mu, const_alpha_s);
	initialize_transport_coeff(A, B);

	AllProcesses[4] = std::vector<Process>();
	AllProcesses[4].push_back( Rate22("Boltzmann/cq2cq", path, dX_Qq2Qq_dt) ); // 2->2
	AllProcesses[4].push_back( Rate22("Boltzmann/cg2cg", path, dX_Qg2Qg_dt) ); // 2->2
	AllProcesses[4].push_back( Rate23("Boltzmann/cq2cqg", path, M2_Qq2Qqg_0) ); // 2->3
	AllProcesses[4].push_back( Rate23("Boltzmann/cg2cgg", path, M2_Qg2Qgg_0) ); // 2->3
	AllProcesses[4].push_back( Rate32("Boltzmann/cqg2cq", path, Ker_Qqg2Qq) );  // 3->2
	AllProcesses[4].push_back( Rate32("Boltzmann/cgg2cg", path, Ker_Qgg2Qg) );  // 3->2
    AllProcesses[4].push_back( Rate12("Boltzmann/c2cg", path, LGV_Q2Qg) );  // 1->2
    AllProcesses[4].push_back( Rate21("Boltzmann/cg2c", path, LGV_Qg2Q) );  // 2->1

    AllProcesses[5] = std::vector<Process>();
    AllProcesses[5].push_back( Rate22("Boltzmann/bq2bq", path, dX_Qq2Qq_dt) );// 2->2
    AllProcesses[5].push_back( Rate22("Boltzmann/bg2bg", path, dX_Qg2Qg_dt) );// 2->2
    AllProcesses[5].push_back( Rate23("Boltzmann/bq2bqg", path, M2_Qq2Qqg_0) );// 2->3
    AllProcesses[5].push_back( Rate23("Boltzmann/bg2bgg", path, M2_Qg2Qgg_0) );// 2->3
	AllProcesses[5].push_back( Rate32("Boltzmann/bqg2bq", path, Ker_Qqg2Qq) );// 3->2
	AllProcesses[5].push_back( Rate32("Boltzmann/bgg2bg", path, Ker_Qgg2Qg) );// 3->2
    AllProcesses[5].push_back( Rate12("Boltzmann/b2bg", path, LGV_Q2Qg) );// 1->2
    AllProcesses[5].push_back( Rate21("Boltzmann/bg2b", path, LGV_Qg2Q) );// 2->1

    AllProcesses[21] = std::vector<Process>();
    AllProcesses[21].push_back( Rate22("Boltzmann/gq2gq", path, dX_gq2gq_dt) ); // 2->2
    AllProcesses[21].push_back( Rate22("Boltzmann/gg2gg", path, dX_gg2gg_dt) ); // 2->2


	BOOST_FOREACH(Process& r, AllProcesses[4]) init_process(r, mode);
	BOOST_FOREACH(Process& r, AllProcesses[5]) init_process(r, mode);
	BOOST_FOREACH(Process& r, AllProcesses[21]) init_process(r, mode);
}

int update_particle_momentum_HT(double dt, double temp, std::vector<double> v3cell, particle & pIn){
	int absid = std::abs(pIn.pid);
	// A Langevin step, does nothing if A and B very small
	fourvec pnew;
	Ito_update(pIn.pid, dt, pIn.mass, temp, v3cell, pIn.p, pnew);
	pIn.p = pnew;
	
	// Matrix-element scattering
	std::vector<fourvec> FS;
	int channel=0.;
	auto p_cell = pIn.p.boost_to(v3cell[0], v3cell[1], v3cell[2]);
	double D_formation_t23_cell = (pIn.x.t()-pIn.t_rad) / pIn.p.t() * p_cell.t();
	double D_formation_t32_cell = (pIn.x.t()-pIn.t_abs) / pIn.p.t() * p_cell.t();
	double dt_cell = dt / pIn.p.t() * p_cell.t();
	double E_cell = p_cell.t();
    std::vector<double> P_channels(AllProcesses[absid].size());
	double P_total = 0., dP;
	double qhatg = qhat_pQCD(21, E_cell, temp);
	BOOST_FOREACH(Process& r, AllProcesses[absid]){
		switch(r.which()){
			case 0:
				if (boost::get<Rate22>(r).IsActive())
					dP = boost::get<Rate22>(r).GetZeroM({E_cell, temp}).s * dt_cell;
				else dP = 0.0;
				P_channels[channel] = P_total + dP;
				break;
			case 1:
				if (boost::get<Rate23>(r).IsActive())
					dP = boost::get<Rate23>(r).GetZeroM({E_cell, temp, D_formation_t23_cell}).s * dt_cell;
				else dP = 0.0;
				P_channels[channel] = P_total + dP;
				break;
			case 2:
				if (boost::get<Rate32>(r).IsActive())
					dP = boost::get<Rate32>(r).GetZeroM({E_cell, temp, D_formation_t32_cell}).s * dt_cell;
				else dP = 0.0;
				P_channels[channel] = P_total + dP;
				break;
			case 3:
				if (boost::get<Rate12>(r).IsActive())
					dP = qhatg*boost::get<Rate12>(r).GetZeroM({E_cell, temp, D_formation_t23_cell}).s * dt_cell;
				else dP = 0.0;
				P_channels[channel] = P_total + dP;
				break;
			case 4:
				if (boost::get<Rate21>(r).IsActive())
					dP = qhatg*boost::get<Rate21>(r).GetZeroM({E_cell, temp, D_formation_t32_cell}).s * dt_cell;
				else dP  = 0.0;
				P_channels[channel] = P_total + dP;
				break;
			break;
				default:
					exit(-1);
					break;
		}
		P_total += dP;
		channel ++;
	}
	for(auto& item : P_channels) {item /= P_total;}
	if (P_total > 0.20) LOG_WARNING << "P_total = " << P_total << " may be too large";
	if ( Srandom::init_dis(Srandom::gen) > P_total) return -1;
	else{
		double p = Srandom::init_dis(Srandom::gen);
		for(int i=0; i<P_channels.size(); ++i){
			if (P_channels[i] > p) {
				channel = i;
				break;
			}
		}
	}
	// Do scattering
	switch(AllProcesses[absid][channel].which()){
		case 0:
			boost::get<Rate22>(AllProcesses[absid][channel]).sample({E_cell, temp}, FS);
			break;
		case 1:
			boost::get<Rate23>(AllProcesses[absid][channel]).sample({E_cell, temp, D_formation_t23_cell}, FS);
			break;
		case 2:
			boost::get<Rate32>(AllProcesses[absid][channel]).sample({E_cell, temp, D_formation_t32_cell}, FS);
			break;
		case 3:
		    boost::get<Rate12>(AllProcesses[absid][channel]).sample({E_cell, temp, D_formation_t23_cell}, FS);
	    	break;
		case 4:
	    	boost::get<Rate21>(AllProcesses[absid][channel]).sample({E_cell, temp, D_formation_t32_cell}, FS);
	    	break;
		default:
			LOG_FATAL << "Channel = " << channel << " not exists";
			exit(-1);
			break;
	}
	// rotate it back and boost it back
	for(auto & pmu : FS) {
		pmu = pmu.rotate_back(p_cell);
		pmu = pmu.boost_back(v3cell[0], v3cell[1], v3cell[2]);
	}
	// Update formation time
	if (channel == 2 || channel == 3 || channel == 6) pIn.t_rad = pIn.x.t();
	if (channel == 4 || channel == 5 || channel == 7) pIn.t_abs = pIn.x.t();

	return channel;
}

// Gluon elastic scattering
int gluon_elastic_scattering(double dt, double temp, std::vector<double> v3cell, fourvec incoming_p, fourvec & outgoing_p){
	int absid = 21;
	// A Langevin step, does nothing if A and B very small
	fourvec pnew;
	Ito_update(absid, dt, 0.0, temp, v3cell, incoming_p, pnew);
	incoming_p = pnew;
	
	// Elastic Scattering
	std::vector<fourvec> FS;
	int channel = 0;
	auto p_cell = incoming_p.boost_to(v3cell[0], v3cell[1], v3cell[2]);
	double dt_cell = dt / incoming_p.t() * p_cell.t();
	double E_cell = p_cell.t();
    std::vector<double> P_channels(AllProcesses[absid].size());
	double P_total = 0., dP;
	BOOST_FOREACH(Process& r, AllProcesses[absid]){
		switch(r.which()){
			case 0:
				if (boost::get<Rate22>(r).IsActive())
					dP = boost::get<Rate22>(r).GetZeroM({E_cell, temp}).s * dt_cell;
				else dP = 0.0;
				P_channels[channel] = P_total + dP;
				break;
			default:
				exit(-1);
				break;
		}
		P_total += dP;
		channel ++;
	}
	for(auto& item : P_channels) {item /= P_total;}
	if (P_total > 0.15 && !type1_warned) {
		LOG_WARNING << "P(g) = " << P_total << " may be too large";
		type1_warned = true;
	}
	if ( Srandom::init_dis(Srandom::gen) > P_total) {
		outgoing_p = incoming_p;
		return -1;
	}
	else{
		double p = Srandom::init_dis(Srandom::gen);
		for(int i=0; i<P_channels.size(); ++i){
			if (P_channels[i] > p) {
				channel = i;
				break;
			}
		}
	}
	// Do scattering
	switch(AllProcesses[absid][channel].which()){
		case 0:
			boost::get<Rate22>(AllProcesses[absid][channel]).sample({E_cell, temp}, FS);
			break;
		default:
			LOG_FATAL << "g-Channel " << channel << " does not exist";
			exit(-1);
			break;
	}
	// rotate it back and boost it back
	for(auto & pmu : FS) {
		pmu = pmu.rotate_back(p_cell);
		pmu = pmu.boost_back(v3cell[0], v3cell[1], v3cell[2]);
	}

	outgoing_p = FS[0];
	return channel;
}

double formation_time(fourvec p, fourvec k, double M, double T, double mpf){
	double x = k.t()/p.t();
	double mD2 = t_channel_mD2->get_mD2(T);
	double mg2 = mD2/2.;
	double Q0 = 2*x*(1-x)*p.t();
	double mass_sqrs = x*x*M*M + (1.-x)*mg2;

	double QT2_LL = measure_perp(p,k).pabs2()*(1.-x+4./9.*x*x);
	double tauf_LL = Q0/(QT2_LL + mass_sqrs);
	return tauf_LL;
}

std::ofstream fc("stat.dat");

int update_particle_momentum_Lido(double dt, double temp, std::vector<double> v3cell, particle & pIn){
	int absid = std::abs(pIn.pid);
	// A Langevin step, does nothing if A and B very small
	fourvec pnew;
	Ito_update(pIn.pid, dt, pIn.mass, temp, v3cell, pIn.p, pnew);
	pIn.p = pnew;

	// Maxtrix element scattering
	std::vector<fourvec> FS; 
	int channel = 0;
	auto p_cell = pIn.p.boost_to(v3cell[0], v3cell[1], v3cell[2]);
	double dt_cell = dt / pIn.p.t() * p_cell.t();
	double E_cell = p_cell.t();
    std::vector<double> P_channels(AllProcesses[absid].size());
	double P_total = 0., dP;

	// Calculate Rate for Heavy Quark
	double qhatg = qhat_pQCD(21, E_cell, temp);
	BOOST_FOREACH(Process& r, AllProcesses[absid]){
		switch(r.which()){
			case 0:
				if (boost::get<Rate22>(r).IsActive()){
					dP = boost::get<Rate22>(r).GetZeroM({E_cell, temp}).s * dt_cell;
				}
				else {dP = 0.0;}
				P_channels[channel] = P_total + dP;
				break;
			case 1:
				if (boost::get<Rate23>(r).IsActive()){
					dP = boost::get<Rate23>(r).GetZeroM({E_cell, temp}).s * dt_cell;
				}
				else {dP = 0.0;}
				P_channels[channel] = P_total + dP;
				break;
			case 2:
				if (boost::get<Rate32>(r).IsActive())
					dP = boost::get<Rate32>(r).GetZeroM({E_cell, temp}).s * dt_cell;
				else dP = 0.0;
				P_channels[channel] = P_total + dP;
				break;
			case 3:
				if (boost::get<Rate12>(r).IsActive())
					dP = qhatg*boost::get<Rate12>(r).GetZeroM({E_cell, temp}).s * dt_cell;
				else dP = 0.0;
				P_channels[channel] = P_total + dP;
				break;
			case 4:
				if (boost::get<Rate21>(r).IsActive())
					dP = qhatg*boost::get<Rate21>(r).GetZeroM({E_cell, temp}).s * dt_cell;
				else dP  = 0.0;
				P_channels[channel] = P_total + dP;
				break;
			default:
				LOG_FATAL << "1. ProcessType = " << r.which() << " not exists";
				exit(-1);
				break;
		}
		P_total += dP;
		channel ++;
	}
	// Normalize Rate
	
	for(auto& item : P_channels) {item /= P_total;}
	if (P_total > 0.15 && !type2_warned) {
		LOG_WARNING << "P(Q) = " << P_total << " may be too large";
		type2_warned = true;
	}
	// Sample channel
	if ( Srandom::init_dis(Srandom::gen) > P_total) channel = -1;
	else{
		double p = Srandom::init_dis(Srandom::gen);
		for(int i=0; i<P_channels.size(); ++i){
			if (P_channels[i] > p) {
				channel = i;
				break;
			}
		}
	}

	if (channel >= 0){
		// Do scattering
		switch(AllProcesses[absid][channel].which()){
			case 0:
				boost::get<Rate22>(AllProcesses[absid][channel]).sample({E_cell, temp}, FS);
				break;
			case 1:
				boost::get<Rate23>(AllProcesses[absid][channel]).sample({E_cell, temp}, FS);
				break;
			case 2:
				boost::get<Rate32>(AllProcesses[absid][channel]).sample({E_cell, temp}, FS);
				break;
			case 3:
				boost::get<Rate12>(AllProcesses[absid][channel]).sample({E_cell, temp}, FS);
				break;
			case 4:
				boost::get<Rate21>(AllProcesses[absid][channel]).sample({E_cell, temp}, FS);
				break;
			default:
				LOG_FATAL << "2. Channel = " << channel << " not exists";
				exit(-1);
				break;
		}
		// rotate it back and boost it back to lab frame
		for(auto & pmu : FS) {
			pmu = pmu.rotate_back(p_cell);
			pmu = pmu.boost_back(v3cell[0], v3cell[1], v3cell[2]);
		}
		// elastic processes change it inmediately
		if (channel ==0 || channel == 1) pIn.p = FS[0]; 

		if(channel == 2 || channel == 3){
			// If it radiates, add a pre-gluon
			pregluon g;
			g.p0 = pIn.p;
			g.k1 = FS[2]; 
			g.kn = g.k1;
			g.t0 = pIn.x.t();
			g.T0 = temp;
			
			double xfrac = g.k1.t()/g.p0.t();
			double local_qhat = 0.;
			double k1_cell = g.k1.boost_to(v3cell[0], v3cell[1], v3cell[2]).t();
			BOOST_FOREACH(Process& r, AllProcesses[21]){
				switch(r.which()){
					case 0:
						if (boost::get<Rate22>(r).IsActive()){
							tensor A = boost::get<Rate22>(r).GetSecondM(
								{std::min(k1_cell, E_cell-k1_cell), temp}
							);			
							local_qhat += A.T[1][1]+A.T[2][2];
						}
						break;
					default:
						break;
				}
			}
			double mD2 = t_channel_mD2->get_mD2(temp);
			g.local_mfp = 1.5*mD2/local_qhat*g.k1.t()/k1_cell; // transform back to lab frame
			pIn.radlist.push_back(g);
		}
		if(channel == 6){
			// If it radiates, add a pre-gluon
			pregluon g;
			g.p0 = pIn.p;
			g.k1 = FS[1]; 
			g.kn = g.k1;
			g.t0 = pIn.x.t();
			g.T0 = temp;
			
			double local_qhat = qhat_pQCD(21, E_cell, temp);
			double mD2 = t_channel_mD2->get_mD2(temp);
			g.local_mfp = .5*mD2/local_qhat*pIn.p.t()/E_cell; // transform back to lab frame
			pIn.radlist.push_back(g);
		}
		if (channel == 4 || channel == 5){
			// If it absorb, add a pre-gluon
			pregluon g;
			g.p0 = pIn.p;
			g.k1 = FS[2]; 
			g.kn = g.k1;
			g.t0 = pIn.x.t();
			g.T0 = temp;
			
			double xfrac = g.k1.t()/g.p0.t();
			double local_qhat = 0.;
			double k1_cell = g.k1.boost_to(v3cell[0], v3cell[1], v3cell[2]).t();
			BOOST_FOREACH(Process& r, AllProcesses[21]){
				switch(r.which()){
					case 0:
						if (boost::get<Rate22>(r).IsActive()){
							tensor A = boost::get<Rate22>(r).GetSecondM(
								{std::min(k1_cell, (1-xfrac)*E_cell), temp}
							);			
							local_qhat += A.T[1][1]+A.T[2][2];
						}
						break;
					default:
						break;
				}
			}
			double mD2 = t_channel_mD2->get_mD2(temp);
			g.local_mfp = mD2/local_qhat*g.k1.t()/k1_cell; // transform back to lab frame
			pIn.abslist.push_back(g);
		}
		if(channel == 7){
			// If it absorb, add a pre-gluon
			pregluon g;
			g.p0 = pIn.p;
			g.k1 = FS[1]; 
			g.kn = g.k1;
			g.t0 = pIn.x.t();
			g.T0 = temp;
			// add diffusion part
			double local_qhat = qhat_pQCD(21, E_cell, temp);
			double mD2 = t_channel_mD2->get_mD2(temp);
			g.local_mfp = mD2/local_qhat*pIn.p.t()/E_cell; // transform back to lab frame
			pIn.abslist.push_back(g);
		}
		if(channel>= AllProcesses[absid].size()){
			LOG_FATAL << "3. Channel = " << channel << " not exists";
			exit(-1);
		}
	}

	// Perform pre-gluon scattering if there is one
	if (!pIn.radlist.empty()){
		// loop over each pre-gluon
		for(std::vector<pregluon>::iterator it=pIn.radlist.begin(); it!=pIn.radlist.end();){
			double taun = formation_time(it->p0,it->kn,pIn.mass,temp,it->local_mfp);
			if (pIn.x.t()-it->t0 > taun){
				double kt20 = measure_perp(it->p0, it->k1).pabs2();
				double kt2n = measure_perp(it->p0, it->kn).pabs2();
				double theta2 = kt2n/std::pow(it->kn.t(),2);
				double thetaM2 = std::pow(pIn.mass/it->p0.t(),2);
				double mD2 = t_channel_mD2->get_mD2(temp);

				double LPM = it->local_mfp/taun
					*std::sqrt(  std::log(1+taun/it->local_mfp)
								/std::log(1+6*it->k1.t()*temp/mD2) );
				double DeadCone = std::pow(theta2/(theta2+thetaM2), 4);
				double RuningCoupling = alpha_s(kt2n, it->T0)/alpha_s(kt20, it->T0);
				double Acceptance = std::min(1.0, LPM*RuningCoupling*DeadCone);


				if (Srandom::rejection(Srandom::gen) < Acceptance){
					pIn.p = pIn.p - it->k1;
					pIn.p.a[0] = std::sqrt(pIn.p.pabs2()+pIn.mass*pIn.mass);
					fc << it->t0 << " " << it->k1.t() << " " << kt20 << " " << kt2n << " " << taun << std::endl;
					
				}
				it = pIn.radlist.erase(it);
			}else{ // else, evolve it
				fourvec k_new;
				gluon_elastic_scattering(dt, temp, v3cell, it->kn, k_new);
				it->kn = k_new*(it->kn.t()/k_new.t());
				it++;
			}
		}
	}

	if (!pIn.abslist.empty()){
		// loop over each pre-gluon
		for(std::vector<pregluon>::iterator it=pIn.abslist.begin(); it!=pIn.abslist.end();){
			fourvec ptot = it->p0+it->kn;
			ptot.a[0] = std::sqrt(ptot.pabs2()+pIn.mass*pIn.mass);
			double taun = formation_time(ptot,it->kn,pIn.mass,temp,it->local_mfp);
			if (pIn.x.t()-it->t0 > taun){
				double kt20 = measure_perp(it->p0, it->k1).pabs2();
				double kt2n = measure_perp(it->p0, it->kn).pabs2();
				double theta2 = kt2n/std::pow(it->kn.t(),2);
				double thetaM2 = std::pow(pIn.mass/(it->p0.t()+it->kn.t()),2);
				double mD2 = t_channel_mD2->get_mD2(temp);

				double LPM = it->local_mfp/taun
					*std::sqrt(  std::log(1+taun/it->local_mfp)
								/std::log(1+6*it->k1.t()*temp/mD2) );
				double DeadCone = std::pow(theta2/(theta2+thetaM2), 4);
				double RuningCoupling = alpha_s(kt2n, it->T0)/alpha_s(kt20, it->T0);
				double Acceptance = std::min(1.0, LPM*RuningCoupling*DeadCone);
				if (Srandom::rejection(Srandom::gen) < Acceptance){
					pIn.p = pIn.p + it->k1;
					pIn.p.a[0] = std::sqrt(pIn.p.pabs2()+pIn.mass*pIn.mass);
				}
				it = pIn.abslist.erase(it);
			}else{ // else, evolve it
				fourvec k_new;
				gluon_elastic_scattering(dt, temp, v3cell, it->kn, k_new);
				it->kn = k_new*(it->kn.t()/k_new.t());
				it++;
			}
		}
	}

	return channel;
}

std::vector<double> probe_test(double M, double E0, double T, double dt, int Nsteps, int Nparticles, std::string mode, double mu, double const_alphas, double A, double B){
	initialize(mode, "./settings.xml", mu, const_alphas, A, B);
	double fmc_to_GeV_m1 = 5.026;
	std::vector<double> dE;
	std::vector<particle> plist(Nparticles);
	double pabs0 = std::sqrt(E0*E0-M*M);
	fourvec p0{E0, 0, 0, pabs0};
	for (auto & p : plist) {
		p.mass = M;
		p.pid = 4;
		p.x = fourvec{0,0,0,0};
		p.p = p0;
		p.t_rad = 0.;
		p.t_abs = 0.;
	}
	double time = 0.;
    double sum = 0.;
	for (int it=0; it<Nsteps; ++it){
		if (it%100 ==0) LOG_INFO << it << " steps, " << "time = " << time << " [fm/c]";
		for (auto & p : plist) sum += E0-p.p.t();
		dE.push_back(sum/Nparticles);
		time += dt;

		for (auto & p : plist){
      		// reset energy for probe test
			p.p = p0;
			// x-update
			p.freestream(dt*fmc_to_GeV_m1);
			// p-update
			int channel = update_particle_momentum_Lido(dt*fmc_to_GeV_m1, T, {0.0, 0.0, 0.0}, p);
		}

	}
	return dE;
}

std::vector<double> Bjorken_test(double M, double E0, double T0, double t0, double dt, int Nsteps, int Nparticles, std::string mode, double mu, double const_alphas, double A, double B){
	initialize(mode, "./settings.xml", mu, const_alphas, A, B);
	double fmc_to_GeV_m1 = 5.026;
	std::vector<double> dE;
	std::vector<particle> plist(Nparticles);
	double pabs0 = std::sqrt(E0*E0-M*M);
	fourvec p0{E0, 0, 0, pabs0};

	for (auto & p : plist) {
		p.pid = 4;
		p.x = fourvec{0,0,0,0};
		p.p = p0;
		p.t_rad = 0.;
		p.t_abs = 0.;
		p.mass = M;
	}
	double time = 0.;
    double sum = 0.;
	for (int it=0; it<Nsteps; ++it){
		double T = T0*std::pow(t0/(t0+time), 1./3.);
		if (it%100 ==0) LOG_INFO << it << " steps, " << "time = " << time << " [fm/c]";
		for (auto & p : plist) sum += E0-p.p.t();
		dE.push_back(sum/Nparticles);
		time += dt;
		for (auto & p : plist){
      		// reset energy for probe test
			p.p = p0;
			// x-update
			p.freestream(dt*fmc_to_GeV_m1);
			// p-update
			int channel = update_particle_momentum_Lido(dt*fmc_to_GeV_m1, T, {0.0, 0.0, 0.0}, p);
		}
	}
	return dE;
}

