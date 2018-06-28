#include <string>
#include <iostream>
#include <exception>

#include "simpleLogger.h"
#include "workflow.h"


void test_config(void);
void test_table(void);

// This sample program evolve heavy quark (E0=30, M=1.3) in a static medium for t=3.0fm/c and calculate the average energy loss (<E-E0>)
int main(int argc, char* argv[]){
	double fmc_to_GeV_m1 = 5.026, M=1.3, E0=30., T=0.3;
	double time = 0.; // system time start from t=0
    double dt = 0.02; // time step 0.02 fm/c
	int Nsteps=150, Nparticles=10000;
	if (argc==1){
		std::cout << "Please tell the program whether use old table (if exists)." << std::endl;
		std::cout << "   $>./example1 new" << std::endl;
		std::cout << "Or $>./example1 old" << std::endl;
		return 1;
	}
	std::string mode = argv[1];
	initialize(mode, "./settings.xml", 1.0);
	std::vector<particle> plist(Nparticles); // a list of particles
	std::vector<fourvec> FS; // final state holder for each scattering
    fourvec p0{E0, 0, 0, std::sqrt(E0*E0-M*M)};
	for (auto & p : plist) {
		p.pid = 4; // pid for charm
		p.x = fourvec{0,0,0,0}; // initialize position at the origin
		p.p = p0; // initial momentum
		p.t_rad = 0.; // initilize the time of last radiation
		p.t_absorb = 0.; // initilize the time of last absorption
	}
	for (int it=0; it<Nsteps; ++it){
		LOG_INFO << it << " steps, " << "time = " << time << " [fm/c]";
		time += dt;
		for (auto & p : plist){
			// loop over each particle
			// perform freestreaming
			p.freestream(dt);
			// Determine scattering channel:
			// <0: nothing happened
			// 0: Qq->Qq
			// 1: Qg->Qg
			// 2: Qq->Qqg
			// 3: Qg->Qgg
			// 4: Qqg->Qq
			// 5: Qgg->Qg
			int channel = update_particle_momentum(
				dt*fmc_to_GeV_m1,  // Delta_t in units of GeV^-1
				T, {0.0, 0.0, 0.0}, // temperature and meidum {vx, vy, vz}
				p.pid, // particle id
				(p.x.t()-p.t_rad)*fmc_to_GeV_m1, // time from last rad [GeV^-1] 
				(p.x.t()-p.t_absorb)*fmc_to_GeV_m1, // time from last rad [GeV^-1] 
				p.p, // particle momentum 
				FS // final state holder
			);
			// Channel 0,1,4,5
			// FS[0] is heavy quark, FS[1] is medium recoil particle
			// Channel 2,3: 
			// FS[0] is heavy quark, FS[1] is medium recoil particle,
			// FS[2] is radiation gluon
			 		 
			// if a scattering happens, update particle momentum and update
			if (channel>=0) {
				p.p = FS[0];
				if (channel == 2 || channel ==3) p.t_rad = time;
				if (channel == 4 || channel ==5) p.t_absorb = time;
			}
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
