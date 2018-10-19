#include <string>
#include <iostream>
#include <exception>

#include "Pythia8/Pythia.h"
#include "simpleLogger.h"
#include "workflow.h"

using namespace Pythia8;

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
    ////////////////////////////////////////////
	double mu = 2.0; 
	double const_alphas = 0.3;
	double A = 1. - 0.25/std::log(2); 
	double B = 0.0; 
	////////////////////////////////////////////
	double M = 1.5; // GeV
	double T = 0.3; // GeV
	double L = 3; // fm
	double dt = 0.02; // fm/c
	double fmc_to_GeV_m1 = 5.026;

	for (int iEvent = 0; iEvent < Neve; ++iEvent) 
    {
	    pythia.next();
		double weight = pythia.info.weight();
        std::map<int,std::vector<std::vector<double>>> gluon_list;
		for (size_t i = 0; i < pythia.event.size(); ++i) 
		{
			auto & pi = pythia.event[i];
     		if(pi.id() == 21 
				&& std::abs(pi.status()) == 51
				)
            {
            	int mother1=pi.mother1();
           		int mother2=pi.mother2();
				auto & m1 = pythia.event[mother1];
            	if(mother2 == 0 && (abs(m1.id())==4))
            	{
                	int motherid=m1.iTopCopyId();
          			gluon_list[motherid].push_back({pi.e(),pi.px(),pi.py(),pi.pz(),
										  			m1.e(),m1.px(),m1.py(),m1.pz()});
            	}
            }
    	}

		for(auto& kv : gluon_list)
		{
		  particle p;
		  auto & Qi = pythia.event[kv.first];
		  double e0=sqrt(Qi.px()*Qi.px()+Qi.py()*Qi.py()+Qi.pz()*Qi.pz()+M*M);
		  fourvec p0{e0, Qi.px(), Qi.py(), Qi.pz()};
		  p.mass = M; // mass
		  p.pid = 4; // charm quark
		  p.x = fourvec{0,0,0,0}; // initial position
		  p.p = p0; // initial momenta
		  p.weight = weight;
		  for (auto& p4 : kv.second)
		  {
		    pregluon G1;
			double qe0=sqrt(p4[5]*p4[5]+p4[6]*p4[6]+p4[7]*p4[7]+M*M);
			G1.p0 = fourvec{qe0,p4[5],p4[6],p4[7]}; // the mother parton energy
			G1.k1 = fourvec{p4[0],p4[1],p4[2],p4[3]}; // the initial four-momenta of the preformed gluon
			G1.kn = G1.k1; // the current four-momenta of the preformed gluon
			G1.t0 = p.x.t(); // creation time, in units of [GeV^{-1}]
			G1.T0 = T; // Temperature at the creation time;
			//G1.local_mfp = ...; // for medium-incuded gluon, this shoul be ~ qhat/mD^2
						// for vacuum radiated gluon, this quantity should not be used
			//Add it to the preformed gluon list
			G1.is_vac=true;
			p.radlist.push_back(G1);
		  }
		  plist.push_back(p);
		}
  }
	
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
	LOG_INFO << "Final pT energy: " << Ef << " GeV";

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


