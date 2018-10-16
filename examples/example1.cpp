#include <string>
#include <iostream>
#include <exception>

#include "simpleLogger.h"
#include "workflow.h"

// This sample program evolve heavy quark (E0=30, M=1.3) in a static medium for t=3.0fm/c and calculate the average energy loss (<E-E0>)
int main(int argc, char* argv[]){
	////////////////////////////////////////////
	double mu = 2.0; // don't change...
	double const_alphas = 0.3; // don't change...
	double A = 1.0; // don't change...
	double B = 0.0; // don't change...
	////////////////////////////////////////////
	double M = 1.5; // GeV
	double T = 0.3; // GeV
	double L = 3; // fm
	double dt = 0.05; // fm/c
	double fmc_to_GeV_m1 = 5.026;
	int Nsteps = int(L/dt);
	

	initialize("old", "./settings.xml", mu, const_alphas, A, B);

	// Initialization
	double E0 = 100; // GeV
	int Nparticles = 10000;
	std::vector<particle> plist(Nparticles);
	double pabs0 = std::sqrt(E0*E0-M*M);
	fourvec p0{E0, 0, 0, pabs0};
	for (auto & p : plist) {
		p.mass = M; // mass
		p.pid = 4; // charm quark
		p.x = fourvec{0,0,0,0}; // initial position
		p.p = p0; // initial momenta
		// uncomment this and modify to initialize the pre-gluon list 
		// create a preformed gluon
		pregluon G1;
		G1.is_vac = true; // 0 for vacuum-shower gluon and 1 for medium-induced gluon
		G1.p0 = p0; // the mother parton energy
		G1.k1 = fourvec{15,0,1,std::sqrt(15*15-1*1)} ; // the initial four-momenta of the preformed gluon
		G1.kn = G1.k1; // the current four-momenta of the preformed gluon
		G1.t0 = p.x.t(); // creation time, in units of [GeV^{-1}]
		G1.T0 = T; // Temperature at the creation time;
		// "a tricky mean-free-path"
		// 1) For medium-induced gluon, this is ~ mD2/qhat
		// 	  and do suppression by rejection on mfp/tau_f_self_consistent
		// 2) For vacuum-shower gluon, a reasonble choice is tau_f_0
		//	  and avoid double counting the phase-space by rejection on
		// 	  tau_f_self_consistent/tau_f_0.
		G1.local_mfp = 0.0;//
		// Add it to the preformed gluon list
		p.radlist.push_back(G1);
		// You can add more preformed gluons
	}
	double time = 0.;
    double sum = 0.;
	
	for (int it=0; it<Nsteps; ++it){
		if (it%10 ==0) LOG_INFO << it << " steps, " << "time = " << time << " [fm/c]";
		time += dt;
		for (auto & p : plist){
			// x-update
			p.freestream(dt*fmc_to_GeV_m1);
			// p-update
			int channel = update_particle_momentum_Lido(dt*fmc_to_GeV_m1, T, {0.0, 0.0, 0.0}, p);
		}
	}
	
	double Ef = 0.;
	for (auto & p : plist) Ef += p.p.t();
	Ef /= Nparticles;

	LOG_INFO << "Initial energy: " << E0 << " GeV";
	LOG_INFO << "Final average energy: " << Ef << " GeV";
	return 0;
}


