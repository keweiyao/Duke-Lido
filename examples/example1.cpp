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
	double A = std::log(2)-0.25; // don't change...
	double B = 0.0; // don't change...
	////////////////////////////////////////////
	double M = 1.3; // GeV
	double T = 0.3; // GeV
	double L = 3; // fm
	double dt = 0.01; // fm/c
	double fmc_to_GeV_m1 = 5.026;
	int Nsteps = int(L/dt);
	

	initialize("new", "./settings.xml", mu, const_alphas, A, B);

	// Initialization
	double E0 = 100; // GeV
	int Nparticles = 1000;
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


