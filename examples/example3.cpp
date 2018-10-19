#include <string>
#include <iostream>
#include <exception>

#include "Pythia8/Pythia.h"
#include "simpleLogger.h"
#include "workflow.h"

using namespace Pythia8;

// This sample program evolve heavy quark (E0=30, M=1.3) in a static medium for t=3.0fm/c and calculate the average energy loss (<E-E0>)
int main(int argc, char* argv[]){
	Pythia pythia;
	Event& event = pythia.event;
	// Declare Pythia object
  	pythia.readString("HardQCD:all = off"); // will repeat this line in the xml for demonstration  
	pythia.readString("HardQCD:hardccbar = on");
  	//pythia.readString("SoftQCD:all = on");
  	pythia.readString("HadronLevel:Decay = off");
  	pythia.readString("HadronLevel:all = off");
    pythia.readString("PartonLevel:all = on");
  	pythia.readString("PartonLevel:ISR = on");
  	pythia.readString("PartonLevel:MPI = off");
  	pythia.readString("PartonLevel:FSR = on");
	pythia.readString("PartonLevel:remnants = off");
	//pythia.readString("PromptPhoton:all=on");
  	pythia.readString("WeakSingleBoson:all=off");
  	pythia.readString("WeakDoubleBoson:all=off");
  	pythia.readString("Beams:idA = 2212");
  	pythia.readString("Beams:idB = 2212");
  	pythia.readString("Beams:eCM = 5020");
  	//pythia.readString("HardQCD:hardbbbar = on");
  	//pythia.readString("HardQCD:hardccbar = on");
  	pythia.readString("PhaseSpace:pTHatMin = 200");
    pythia.readString("PhaseSpace:pTHatMin = 210");
	pythia.readString("PhaseSpace:bias2Selection = off");
	pythia.readString("PhaseSpace:bias2SelectionPow = 6.");
	pythia.readString("Random:setSeed = on");
	pythia.readString("Random:seed = 120");
    pythia.init();

	std::vector<particle> plist;
    ////////////////////////////////////////////
	double mu = 2.0; // don't change...
	double const_alphas = 0.3; // don't change...
	double A = 1.0; // don't change...
	double B = 0.0; // don't change...
	////////////////////////////////////////////
	double M = 1.5; // GeV
	double T = 0.3; // GeV
	double L = 3; // fm
	double dt = 0.02; // fm/c
	double fmc_to_GeV_m1 = 5.026;

	for (int iEvent = 0; iEvent < 1000000; ++iEvent) 
    {

	    pythia.next();
		for (size_t i = 0; i < pythia.event.size(); ++i) 
		{
			auto p = event[i];
			if (p.isFinal() && p.idAbs() == 4) 
			{
				// final momenta 
				fourvec p0{p.e(), p.px(), p.py(), p.pz()};
				particle c_entry;
				c_entry.mass = M; // mass
				c_entry.pid = 4; // charm quark
				c_entry.x = fourvec{0,0,0,0}; // initial position
				// trace back to its first production to 
				// find out final state gluon radiation
				while(true)
				{
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
					if (event[p.daughter1()].idAbs() == 21) 
						gluon_eid = p.daughter1();
					else if (event[p.daughter2()].idAbs() == 21) 
						gluon_eid = p.daughter2();

					if (gluon_eid >0 ) { // a radiated gluon:
						// make it pre-formed gluon
						auto g = event[gluon_eid];
						fourvec k{g.e(), g.px(), g.py(), g.pz()};
						p0 = p0 + k; // add back the pre-gluon, put on shell
						p0.a[0] = std::sqrt(p0.pabs2() + M*M);
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
				}
				c_entry.p = p0; // without FSR energy loss
				plist.push_back(c_entry); 
			}
    	}
  }
	
	int Nsteps = int(L/dt);
	int Nparticles=plist.size();
	double Ei= 0.;
	for (auto & p : plist) Ei += p.p.t();
	Ei /= Nparticles;

	initialize("new", "./settings.xml", mu, const_alphas, A, B);

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

	LOG_INFO << "Initial energy: " << Ei << " GeV";
	LOG_INFO << "Final average energy: " << Ef << " GeV";

	/*
	std::fstream fs;
    fs.open ("hq_gluon_modified.dat", std::fstream::out | std::fstream::app);
    for (auto & p : plist)
	{
		for (auto& g : p.radlist)
		{
			if(g.is_vac == true)
			{
				fs<<g.kn.x()<<" "<<g.kn.y()<<" "<<g.kn.z()<<" "<<g.kn.t()<<" "<<g.p0.x()<<" "<<g.p0.y()<<" "<<g.p0.z()<<" "<<g.p0.t()<<endl;
			}
		}
	}
	*/
	return 0;
}


