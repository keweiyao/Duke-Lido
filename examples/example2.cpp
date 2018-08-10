#include <string>
#include <iostream>
#include <exception>

#include "simpleLogger.h"
#include "workflow.h"
// This sample program evolve heavy quark (E0=30, M=1.3) in a Bjorken medium at mid-rapidity for 0.5<t<5.0fm/c and calculate the average energy loss (<E-E0>)
// T = 0.5GeV(0.5/t)^(1/3)
int main(int argc, char* argv[]){
	double fmc_to_GeV_m1 = 5.026, M=1.3, E0=50., T=0.5;	
	double t0 = 0.5;
	double time = t0; // system time start from t=0
    double dt = 0.005; // time step 0.02 fm/c
	int Nsteps=int((5.0-t0)/dt), Nparticles=10000;
	if (argc==1){
		std::cout << "Please tell the program whether use old table (if exists)." << std::endl;
		std::cout << "   $>./example1 new" << std::endl;
		std::cout << "Or $>./example1 old" << std::endl;
		return 1;
	}
	std::string mode = argv[1];
	initialize(mode, "./settings.xml", 1.0, 0.1);
	std::vector<particle> plist(Nparticles); // a list of particles
	std::vector<fourvec> FS; // final state holder for each scattering
    fourvec p0{E0, 0, 0, std::sqrt(E0*E0-M*M)};
	for (auto & p : plist) {
		p.pid = 4; // pid for charm
		p.x = fourvec{t0*fmc_to_GeV_m1,0,0,0}; // initialize position at the origin
		p.p = p0; // initial momentum
		p.has_k_rad = false;
		p.has_k_abs = false;
		p.t_rad = t0*fmc_to_GeV_m1; // initilize the time of last radiation
		p.t_abs = t0*fmc_to_GeV_m1; // initilize the time of last absorption
		p.k_rad = fourvec{0,0,0,0};
		p.k_abs = fourvec{0,0,0,0};
		p.resum_counts = 0;
		p.mass = M;
	}
	for (int it=0; it<Nsteps; ++it){
		LOG_INFO << it << " steps, " << "time = " << time << " [fm/c], T=" << T << " [GeV]";
		time += dt;
		for (auto & p : plist){
			// loop over each particle
			// perform freestreaming
			p.freestream(dt*fmc_to_GeV_m1);
			// Determine scattering channel:
			// <0: nothing happened
			// 0: Qq->Qq
			// 1: Qg->Qg
			// 2: Qq->Qqg
			// 3: Qg->Qgg
			// 4: Qqg->Qq
			// 5: Qgg->Qg
			int channel = update_particle_momentum_BDMPSZ(dt*fmc_to_GeV_m1, T, {0.0, 0.0, 0.0}, p);
			//int channel = update_particle_momentum_HT(dt*fmc_to_GeV_m1, T, {0.0, 0.0, 0.0}, p);
		}
	}
	// Calculate average loss in energy
	double sum=0.;
	for (auto & p : plist){
		sum += p.p.t();
	}
	double E1 = sum/Nparticles;
	// Print summary
	LOG_INFO << "Summary:";
	LOG_INFO << "E0 = " << E0 << " GeV";
	LOG_INFO << "T = " << T << " GeV";
	LOG_INFO << "t-t0 = " << Nsteps*dt << " fm/c";
	LOG_INFO << "Type = charm quark";
	LOG_INFO << "Average E-loss = " << E1-E0 << " GeV";

	return 0;
}
