#include "Pythia8/Pythia.h"
#include <iostream>
#include <fstream>
#include <iomanip> 
#include <memory>
#include <cstdio>
#include <string>
#include <vector>
#include <exception>

#include "simpleLogger.h"
#include "workflow.h"

using namespace Pythia8;

int pythia8_gen(int Nev, std::vector<particle> & plist) {
	plist.clear();

	// Generator. 
	Pythia pythia;
	Event& event = pythia.event;
	Info& info = pythia.info;

	// read settings
	pythia.readFile("pythia-setting-AA.txt");

	// suppress output
	pythia.readString("Init:showProcesses = off");  
    pythia.readString("Init:showMultipartonInteractions = off");  
    pythia.readString("Init:showChangedSettings = off");  
    pythia.readString("Init:showChangedParticleData = off");  
    pythia.readString("Next:numberCount = 1000");  
    pythia.readString("Next:numberShowInfo = 0");  
    pythia.readString("Next:numberShowProcess = 0");  
    pythia.readString("Next:numberShowEvent = 0"); 

	// init pythia
	pythia.init();
	double Mass = 1.5;
	// Generate events.
	auto event_counter = 0;
	while (event_counter < Nev) {
		auto status = pythia.next();
		auto weight = pythia.info.weight();
		event_counter ++;
		for(auto i = 0; i < pythia.event.size(); i ++) {
			auto p = event[i];
			if (p.isFinal() && p.idAbs() == 4 
				//&& std::abs(p.y())<1.0 
				//&& std::abs(p.pT()-100)<10.
			) {
				// final momenta 
				fourvec p0{p.e(), p.px(), p.py(), p.pz()};
				particle c_entry;
				c_entry.mass = Mass; // mass
				c_entry.pid = 4; // charm quark
				c_entry.x = fourvec{0,0,0,0}; // initial position
				// trace back to its first production to 
				// find out final state gluon radiation
				/*while(true){
					auto m1 = p.mother1();
					auto m2 = p.mother2();
					auto d1 = p.daughter1();
					auto d2 = p.daughter2();

					// trace the fermion line back, else break
					if (event[m1].id() == p.id()) 
						p = event[m1];
					else if (event[m2].id() == p.id()) 
						p = event[m2];
					else break;
					
					int gluon_eid = -1;
					if (event[p.daughter1()].idAbs() == 21 && std::abs(event[p.daughter1()].status()) == 51)
						gluon_eid = p.daughter1();
						
					else if (event[p.daughter2()].idAbs() == 21 && std::abs(event[p.daughter2()].status()) == 51) 
						gluon_eid = p.daughter2();

					if (gluon_eid >0 ) { // a radiated gluon:
						// make it pre-formed gluon
						auto g = event[gluon_eid];
						fourvec k{g.e(), g.px(), g.py(), g.pz()};
						p0 = p0 + k; // add back the pre-gluon, put on shell
						p0.a[0] = std::sqrt(p0.pabs2() + Mass*Mass);
						pregluon G1;
						G1.is_vac = true;// vacuum 
						G1.p0 = p0;
						G1.k1 = k;
						G1.kn = G1.k1;
						G1.t0 = 0.0;
						G1.T0 = 0.0; // not used
						G1.local_mfp = 0.0; // not used
						c_entry.radlist.push_back(G1);
					}
				}*/
				c_entry.p = p0; // without FSR energy loss
				c_entry.weight = weight; // event weight
				//double y0 = std::log((p0.t()+p0.z())/(p0.t()-p0.z()));
				//double pT0 = std::sqrt(p0.x()*p0.x() + p0.y()*p0.y());
				//if (std::abs(y0) < 1.){
				plist.push_back(c_entry); 
				//}
			}
		}
	}
	LOG_INFO << info.weightSum() << " " << info.sigmaGen();
	// Done. 
	return 0;
}

// This sample program evolve heavy quark (E0=30, M=1.3) in a static medium for t=3.0fm/c and calculate the average energy loss (<E-E0>)
int main(int argc, char* argv[]){
	////////////////////////////////////////////
	double mu = 2.0; // don't change...
	double const_alphas = 0.3; // don't change...
	double A = 1.-0.25/std::log(2); // don't change...
	double B = 0.0; // don't change...
	////////////////////////////////////////////
	double M = 1.5; // GeV
	double T0 = 0.3; // GeV
	double L = 0.0; // fm
	double dt = .1; // fm/c
	double fmc_to_GeV_m1 = 5.026;
	int Nsteps = int(L/dt);
	

	initialize("old", "./settings.xml", mu, const_alphas, A, B);

	// Initialization
	int Neve = 100000;
	std::vector<particle> plist;
	pythia8_gen(Neve, plist);
	int Ncharm = plist.size();
	
	double Ei = 0.;
	for (auto & p : plist) Ei += std::sqrt(p.p.x()*p.p.x()+p.p.y()*p.p.y());
	Ei /= Ncharm;

	double time = 0.;	
	for (int it=0; it<Nsteps; ++it){
		if (it%10 ==0) LOG_INFO << it << " steps, " << "time = " << time << " [fm/c]";
		time += dt;
		double T = T0;
		if (time > 0.) T = 0.01;
		for (auto & p : plist){
			// x-update
			//auto pi = p.p;
			p.freestream(dt*fmc_to_GeV_m1);
			// p-update
			int channel = update_particle_momentum_Lido(dt*fmc_to_GeV_m1, T, {0.0, 0.0, 0.0}, p);
			//p.p = pi;
		}
	}
	
	double Ef = 0.;
	std::ofstream f("list-AA.dat");
	for (auto & p : plist) {
		Ef += std::sqrt(p.p.x()*p.p.x()+p.p.y()*p.p.y());
		f << p.p << " " << p.weight << std::endl;
	}
	f.close();
	Ef /= Ncharm;

	LOG_INFO << "Initial pT: " << Ei << " GeV";
	LOG_INFO << "Final average pT: " << Ef << " GeV";
	return 0;
}


