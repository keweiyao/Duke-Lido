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
                        default:
                                exit(-1);
                                break;
                }
}


void initialize(std::string mode, std::string path, double mu){
	print_logo();
	boost::property_tree::ptree config;
    initialize_mD_and_scale(1, mu);

	AllProcesses[4] = std::vector<Process>();
	AllProcesses[4].push_back( Rate22("Boltzmann/cq2cq", path, dX_Qq2Qq_dt) );
	AllProcesses[4].push_back( Rate22("Boltzmann/cg2cg", path, dX_Qg2Qg_dt) );
	AllProcesses[4].push_back( Rate23("Boltzmann/cq2cqg", path, M2_Qq2Qqg_0) );
	AllProcesses[4].push_back( Rate23("Boltzmann/cg2cgg", path, M2_Qg2Qgg_0) );
	AllProcesses[4].push_back( Rate32("Boltzmann/cqg2cq", path, Ker_Qqg2Qq) );
	AllProcesses[4].push_back( Rate32("Boltzmann/cgg2cg", path, Ker_Qgg2Qg) );

    AllProcesses[5] = std::vector<Process>();
    AllProcesses[5].push_back( Rate22("Boltzmann/bq2bq", path, dX_Qq2Qq_dt) );
    AllProcesses[5].push_back( Rate22("Boltzmann/bg2bg", path, dX_Qg2Qg_dt) );
    AllProcesses[5].push_back( Rate23("Boltzmann/bq2bqg", path, M2_Qq2Qqg_0) );
    AllProcesses[5].push_back( Rate23("Boltzmann/bg2bgg", path, M2_Qg2Qgg_0) );
	AllProcesses[5].push_back( Rate32("Boltzmann/bqg2bq", path, Ker_Qqg2Qq) );
	AllProcesses[5].push_back( Rate32("Boltzmann/bgg2bg", path, Ker_Qgg2Qg) );

    AllProcesses[21] = std::vector<Process>();
    AllProcesses[21].push_back( Rate22("Boltzmann/gq2gq", path, dX_gq2gq_dt) );
    AllProcesses[21].push_back( Rate22("Boltzmann/gg2gg", path, dX_gg2gg_dt) );


	BOOST_FOREACH(Process& r, AllProcesses[4]) init_process(r, mode);
	BOOST_FOREACH(Process& r, AllProcesses[5]) init_process(r, mode);
	BOOST_FOREACH(Process& r, AllProcesses[21]) init_process(r, mode);
}

int update_particle_momentum_HT(double dt, double temp, std::vector<double> v3cell, particle & pIn){
	std::vector<fourvec> FS;
	int absid = std::abs(pIn.pid), channel=0.;
	auto p_cell = pIn.p.boost_to(v3cell[0], v3cell[1], v3cell[2]);
	double D_formation_t23_cell = (pIn.x.t()-pIn.t_rad) / pIn.p.t() * p_cell.t();
	double D_formation_t32_cell = (pIn.x.t()-pIn.t_abs) / pIn.p.t() * p_cell.t();
	double dt_cell = dt / pIn.p.t() * p_cell.t();
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
			default:
				exit(-1);
				break;
		}
		P_total += dP;
		channel ++;
	}
	for(auto& item : P_channels) {item /= P_total;}
	if (P_total > 0.15) LOG_WARNING << "P_total = " << P_total << " may be too large";
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
	if (channel == 2 || channel == 3) pIn.t_rad = pIn.x.t();
	if (channel == 4 || channel == 5) pIn.t_abs = pIn.x.t();
	pIn.p = FS[0];	
	return channel;
}

// Gluon elastic scattering
int gluon_elastic_scattering(double dt, double temp, std::vector<double> v3cell, fourvec incoming_p, fourvec & outgoing_p){
	std::vector<fourvec> FS;
	int absid = 21, channel = 0;
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
	if (P_total > 0.15) LOG_WARNING << "P(g) = " << P_total << " may be too large";
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

double formation_time(fourvec p, fourvec k, double M, double T){
	// do everything in the frame where heavy quark moves in z-direction
	double p0 = p.t();
	double pz = std::sqrt(p0*p0 - M*M);
	double k0 = k.t();
	double kz = (k.x()*p.x() + k.y()*p.y() + k.z()*p.z() ) / pz;
	double x = (k0+kz)/(p0+pz+k0+kz);
	double kt2 = k0*k0 - kz*kz;
	double mD2 = t_channel_mD2->get_mD2(T);
	double tauf = 2.*x*(1-x)*p0/(kt2 + x*x*M*M + (1-x)*mD2/2.);
	return tauf;
}
//std::ofstream fc("counts.dat");

int update_particle_momentum_BDMPSZ(double dt, double temp, std::vector<double> v3cell, particle & pIn){
	std::vector<fourvec> FS; // Final state holder
	// Perform pre-gluon scattering if there is one
	if (pIn.has_k_rad){
		//LOG_INFO << formation_time(pIn.p, pIn.k_rad, pIn.mass, temp);
		// check if it is formed already (t-t0 > tau_f)
		// if formed, let the heavy quark forget the gluon, else keep evolving it
		if (pIn.x.t()-pIn.t_rad > formation_time(pIn.p, pIn.k_rad, pIn.mass, temp)) {
			pIn.has_k_rad = false;
			//fc << pIn.resum_counts << std::endl;
			// take out the gluon energy from the original heavy quark, note the possible mismatch in 4-momentum conservation due to other rescatterings
			pIn.p.a[1] = pIn.p.x() - pIn.k_rad.x();
			pIn.p.a[2] = pIn.p.y() - pIn.k_rad.y();
			pIn.p.a[3] = pIn.p.z() - pIn.k_rad.z();
			pIn.p.a[0] = std::sqrt(pIn.p.x()*pIn.p.x() + pIn.p.y()*pIn.p.y() + pIn.p.z()*pIn.p.z() + pIn.mass*pIn.mass);
			pIn.resum_counts = 0;
		}	
		else{ 
			fourvec k_rad_new;
			// if the pre-gluon scatterings, update its momentum.
			if ( gluon_elastic_scattering(dt, temp, v3cell, pIn.k_rad, k_rad_new) >= 0) {
				pIn.k_rad = k_rad_new;		
				pIn.resum_counts ++;
			}
		}
	}

	// Now we can preform the heavy quark scattering as always,
	// except that: we don't let the heavy quark radiate if it has a pre-gluon
	int absid = std::abs(pIn.pid), channel = 0;
	auto p_cell = pIn.p.boost_to(v3cell[0], v3cell[1], v3cell[2]);
	double dt_cell = dt / pIn.p.t() * p_cell.t();
	double E_cell = p_cell.t();
    std::vector<double> P_channels(AllProcesses[absid].size());
	double P_total = 0., dP;
	// Calculate Rate for Heavy Quark
	BOOST_FOREACH(Process& r, AllProcesses[absid]){
		switch(r.which()){
			case 0:
				if (boost::get<Rate22>(r).IsActive())
					dP = boost::get<Rate22>(r).GetZeroM({E_cell, temp}).s * dt_cell;
				else dP = 0.0;
				P_channels[channel] = P_total + dP;
				break;
			case 1:
				// if 2->3 is on and there is no pre-gluon, include 2->3 rate
				if (boost::get<Rate23>(r).IsActive() && (!pIn.has_k_rad))
					dP = boost::get<Rate23>(r).GetZeroM({E_cell, temp}).s * dt_cell;
				else dP = 0.0;
				P_channels[channel] = P_total + dP;
				break;
			case 2:
				if (boost::get<Rate32>(r).IsActive())
					dP = boost::get<Rate32>(r).GetZeroM({E_cell, temp}).s * dt_cell;
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
	// Normalize Rate
	for(auto& item : P_channels) {item /= P_total;}
	if (P_total > 0.15) LOG_WARNING << "P_total = " << P_total << " may be too large";
	// Sample Rate
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
			boost::get<Rate23>(AllProcesses[absid][channel]).sample({E_cell, temp}, FS);
			break;
		case 2:
			boost::get<Rate32>(AllProcesses[absid][channel]).sample({E_cell, temp}, FS);
			break;
		default:
			LOG_FATAL << "Channel = " << channel << " not exists";
			exit(-1);
			break;
	}
	// rotate it back and boost it back to lab frame
	for(auto & pmu : FS) {
		pmu = pmu.rotate_back(p_cell);
		pmu = pmu.boost_back(v3cell[0], v3cell[1], v3cell[2]);
	}
	// If it radiates, add a pre-gluon
	if (channel == 2 || channel == 3){
		pIn.t_rad = pIn.x.t();
		pIn.k_rad = FS[2];
		pIn.has_k_rad = true;
		pIn.resum_counts = 1;
		// merge the k into the Q 
		//pIn.p.a[1] = FS[0].x() + pIn.k_rad.x();
		//pIn.p.a[2] = FS[0].y() + pIn.k_rad.y();
		//pIn.p.a[3] = FS[0].z() + pIn.k_rad.z();
		//pIn.p.a[0] = std::sqrt(pIn.p.x()*pIn.p.x() + pIn.p.y()*pIn.p.y() + pIn.p.z()*pIn.p.z() + pIn.mass*pIn.mass);
		//LOG_INFO << pIn.p.t() << " " << pIn.p.x() << " " << pIn.p.y() << " " << pIn.p.z() << " " << dot(pIn.p, pIn.p); 
		//LOG_INFO << pIn.k_rad.t() << " " << pIn.k_rad.x() << " " << pIn.k_rad.y() << " " << pIn.k_rad.z() << " " << dot(pIn.k_rad, pIn.k_rad); 
	}
	else{
		pIn.p = FS[0];
	}
	return channel;
}


std::vector<double> probe_test(double E0, double T, double dt=0.05, int Nsteps=100, int Nparticles=10000, std::string mode="old"){
	double fmc_to_GeV_m1 = 5.026;
	initialize(mode, "./settings.xml", 1.0);
	double M = 1.3;
	std::vector<double> dE;
	std::vector<particle> plist(Nparticles);
    fourvec p0{E0, 0, 0, std::sqrt(E0*E0-M*M)};
	for (auto & p : plist) {
		p.pid = 4;
		p.x = fourvec{0,0,0,0};
		p.p = p0;
		p.has_k_rad = false;
		p.has_k_abs = false;
		p.t_rad = 0.;
		p.t_abs = 0.;
		p.k_rad = fourvec{0,0,0,0};
		p.k_abs = fourvec{0,0,0,0};
		p.mass = M;
	}
	double time = 0.;
    double sum = 0.;
	for (int it=0; it<Nsteps; ++it){
		LOG_INFO << it << " steps, " << "time = " << time << " [fm/c]";
		for (auto & p : plist) sum += E0-p.p.t();
		dE.push_back(sum/Nparticles);
		time += dt;
		for (auto & p : plist){
      		p.p = p0; // reset energy for probe test
			p.freestream(dt*fmc_to_GeV_m1);
			int channel = update_particle_momentum_BDMPSZ(dt*fmc_to_GeV_m1, T, {0.0, 0.0, 0.0}, p);
		}
	}
	return dE;
}

