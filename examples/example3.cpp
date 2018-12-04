#include <string>
#include <iostream>
#include <exception>

#include "Pythia8/Pythia.h"
#include "simpleLogger.h"
#include "workflow.h"

using namespace Pythia8;
std::ofstream f("list-pp-D");
int main(int argc, char* argv[]){
	if (argc < 3){
		LOG_FATAL << "need a pythia setting file and number of event";
		exit(-1);
	}
	std::string filename(argv[1]);
	int Neve = atoi(argv[2]);
	
	// read settings
	Pythia pythia;
	pythia.readFile(filename);
	Event& event = pythia.event;
	Info& info = pythia.info;

	// suppress output
	pythia.readString("SoftQCD:all = off");
	pythia.readString("PromptPhoton:all=off");
  	pythia.readString("WeakSingleBoson:all=off");
  	pythia.readString("WeakDoubleBoson:all=off");

	pythia.readString("Init:showProcesses = off");  
    pythia.readString("Init:showMultipartonInteractions = off");  
    pythia.readString("Init:showChangedSettings = off");  
    pythia.readString("Init:showChangedParticleData = off");  
    pythia.readString("Next:numberCount = 1000");  
    pythia.readString("Next:numberShowInfo = 0");  
    pythia.readString("Next:numberShowProcess = 0");  
    pythia.readString("Next:numberShowEvent = 0"); 

    pythia.init();

	std::vector<particle> plist;
	double M = 1.3; // GeV
	for (int iEvent = 0; iEvent < Neve; ++iEvent) 
    {
	    pythia.next();
		double weight = pythia.info.weight();
		for (size_t i = 0; i < event.size(); ++i) 
		{
			auto p = event[i];
			if (p.isFinal() && 
				(p.idAbs() == 411 || p.idAbs() == 421 || p.idAbs() == 413) )
			{
				f << p.idAbs() << " " << p.e()<< " "<< p.px()<< " "<< p.py()<< " "<< p.pz()<< " "<< weight << std::endl;
				/*// final momenta 
				fourvec p0{p.e(), p.px(), p.py(), p.pz()};
				particle c_entry;
				c_entry.mass = M; // mass
				c_entry.pid = 4; // charm quark
				c_entry.x = fourvec{0,0,0,0}; // initial position
				c_entry.weight = weight;
				// trace back to its first production to 
				// find out final state gluon radiation
				while(true){
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
					if (event[p.daughter1()].idAbs() == 21 
					&& std::abs(event[p.daughter1()].status()) == 51)
						gluon_eid = p.daughter1();
						
					else if (event[p.daughter2()].idAbs() == 21 
					&& std::abs(event[p.daughter2()].status()) == 51) 
						gluon_eid = p.daughter2();

					if (gluon_eid >0) { // a radiated gluon:
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
						G1.t0 = c_entry.x.t();
						G1.T0 = 0.0; // not used
						G1.local_mfp = 0.0; // not used
						c_entry.radlist.push_back(G1);
					}
				}
				c_entry.p = p0; // without FSR energy loss
				plist.push_back(c_entry); */
			}
    	}
  	}
	LOG_INFO << "TW, SG " <<  pythia.info.weightSum() << " " << pythia.info.sigmaGen();
    ////////////////////////////////////////////
	/*double mu = 2.0;
	double const_alphas = 0.3; 
	double A = 1. - 0.25/std::log(2); 
	double B = 0.0; 
	////////////////////////////////////////////
	double T = 0.3; // GeV
	double L = 3; // fm
	double dt = 0.02; // fm/c
	double fmc_to_GeV_m1 = 5.026;
	int Nsteps = int(L/dt);
	int Nparticles=plist.size();


	double Ei= 0.;
	for (auto & p : plist) Ei += std::sqrt(p.p.x()*p.p.x()+p.p.y()*p.p.y());
	Ei /= Nparticles;

	initialize("old", "./settings.xml", mu, const_alphas, A, B);

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
	for (auto & p : plist) Ef += std::sqrt(p.p.x()*p.p.x()+p.p.y()*p.p.y());
	Ef /= Nparticles;

	LOG_INFO << "Nparticles: " << Nparticles;
	LOG_INFO << "Initial pT: " << Ei << " GeV";
	LOG_INFO << "Final pT: " << Ef << " GeV";
	*/

	return 0;
}


