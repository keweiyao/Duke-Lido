#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <exception>
#include <random>

#include "Pythia8/Pythia.h"
#include "simpleLogger.h"
#include "Medium_Reader.h"
#include "workflow.h"


using namespace Pythia8;

class HardGen{
public: 
	HardGen(std::string f_pythia, std::string f_trento, int iev);
	void Generate(std::vector<particle> & plist, int N_pythia_events, 
				int POI, double yabs_cut);
private:
	Pythia pythia;
	const double Mc;
	TransverPositionSampler TRENToSampler;
};

HardGen::HardGen(std::string f_pythia, std::string f_trento, int iev):
Mc(1.3),TRENToSampler(f_trento, iev){	
	// read pythia settings
	pythia.readFile(f_pythia);
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
	// Init
    pythia.init();
	// Read TRENTo Binary collision densisty
}

void HardGen::Generate(std::vector<particle> & plist, int N_pythia_events, 
				int POI, double yabs_cut){
	double x, y;

	for (int iEvent = 0; iEvent < N_pythia_events; ++iEvent) 
    {
	    pythia.next();
		TRENToSampler.SampleXY(x, y);
		double weight = pythia.info.weight();
		for (size_t i = 0; i < pythia.event.size(); ++i) 
		{
			auto p = pythia.event[i];
			if (p.isFinal() && p.idAbs() == POI && std::abs(p.y())<yabs_cut) 
			{
				// final momenta 
				fourvec p0{p.e(), p.px(), p.py(), p.pz()};
				particle c_entry;
				c_entry.mass = Mc; // mass
				c_entry.pid = POI; // charm quark
				c_entry.x = fourvec{0,x,y,0}; // initial position
				c_entry.weight = weight;
				c_entry.freezeout = false;
				c_entry.Tf = 0.0;
				c_entry.vcell.resize(3);
				c_entry.vcell[0] = 0.; 
				c_entry.vcell[1] = 0.; 
				c_entry.vcell[2] = 0.; 
				// trace back to its first production to 
				// find out final state gluon radiation
				while(true){
					auto m1 = p.mother1();
					auto m2 = p.mother2();
					auto d1 = p.daughter1();
					auto d2 = p.daughter2();
					// trace the fermion line back, else break
					if (pythia.event[m1].id() == p.id()) 
						p = pythia.event[m1];
					else if (pythia.event[m2].id() == p.id()) 
						p = pythia.event[m2];
					else break;
					int gluon_eid = -1;
					if (pythia.event[p.daughter1()].idAbs() == 21 
					&& std::abs(pythia.event[p.daughter1()].status()) == 51)
						gluon_eid = p.daughter1();	
					else if (pythia.event[p.daughter2()].idAbs() == 21 
					&& std::abs(pythia.event[p.daughter2()].status()) == 51) 
						gluon_eid = p.daughter2();
					if (gluon_eid >0) { // a radiated gluon:
						// make it pre-formed gluon
						auto g = pythia.event[gluon_eid];
						fourvec k{g.e(), g.px(), g.py(), g.pz()};
						p0 = p0 + k; // add back the pre-gluon, put on shell
						p0.a[0] = std::sqrt(p0.pabs2() + Mc*Mc);
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
				c_entry.p0 = p0;
				plist.push_back(c_entry); 
			}
    	}
  	}

	// get the correct normalization for cross-section
	double normW = pythia.info.sigmaGen()/pythia.info.weightSum();
	for (auto & p : plist) p.weight *= normW;
}

int main(int argc, char* argv[]){
	if (argc < 6){
		std::cout << "Usage: hydro-couple <pythia-setting> <initial-file> <TRTENTo-eid> <hydro-file>  <NpythiaEevents>" << std::endl;
		exit(-1);
	}
	std::string pythia_config(argv[1]);
	std::string trento_initial(argv[2]);
	int trento_eid = atoi(argv[3]);
	std::string hydro_data(argv[4]);
	int NPythiaEvents = atoi(argv[5]);

	std::vector<particle> plist;

	/// HardGen
	HardGen hardgen(pythia_config, trento_initial, trento_eid);
	hardgen.Generate(plist, NPythiaEvents, 4, 2.5);
	
	/// Read Hydro
	Medium<2> med1(hydro_data);
	
	/// Assign each quark a transverse position according to TRENTo Nbin output
	/// Freestream particle to the start of hydro
	for (auto & p : plist){
		double dt_fs = med1.get_tauH()/std::sqrt(1. - p.p.z()*p.p.z()/p.p.t()/p.p.t());
		p.freestream(dt_fs);
	}

    ////////////////////////////////////////////
	double mu = 2.0;
	double const_alphas = 0.3; 
	double A = 1. - 0.25/std::log(2); 
	double B = 0.0; 
	int Nparticles=plist.size();
	initialize("old", "./settings.xml", mu, const_alphas, A, B);

	// initial pT
	double pTi = mean_pT(plist);

	// run
	int counter = 0;
	int Ns = 10;
	std::ofstream f("tab.dat");
	
	while(med1.load_next()){
		double current_hydro_clock = med1.get_tauL();
		double hydro_dtau = med1.get_hydro_time_step();
		for (int i=0; i<Ns; ++i){
			double dtau = hydro_dtau/Ns; // use smaller dt step
			for (auto & p : plist){
				if (p.freezeout) continue; // do not touch freezeout ones
				// determine dt needed to evolve to the next tau
				double tau = std::sqrt(p.x.t()*p.x.t()-p.x.z()*p.x.z());
				double dt_lab = calcualte_dt_from_dtau(p.x, p.p, tau, dtau);
				// get hydro information
				double T = 0., vx = 0., vy = 0., vz = 0.;
				med1.interpolate(p.x, T, vx, vy, vz);
				// x-update
				p.freestream(dt_lab);
				// p-update
				int channel = 
					update_particle_momentum_Lido(dt_lab, T, {vx, vy, vz}, p);
			}
		}
		counter ++;
	}
	f.close();

	// final pT
	double pTf = mean_pT(plist);
	LOG_INFO << "Nparticles: " << Nparticles;
	LOG_INFO << "Initial pT: " << pTi << " GeV";
	LOG_INFO << "Final pT: " << pTf << " GeV";

	output_oscar(plist, "c-quark-frzout.dat");
	return 0;
}



