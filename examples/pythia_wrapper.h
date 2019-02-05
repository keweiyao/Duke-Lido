#ifndef PYTHIA_WRAPPER_H
#define PYTHIA_WRAPPER_H

#include "simpleLogger.h"
#include "Pythia8/Pythia.h"

using namespace Pythia8;

class HardGen{
public: 
    HardGen(std::string f_pythia, std::string f_trento, int iev);
    void Generate(std::vector<particle> & plist, int N_pythia_events, 
                int POI, double yabs_cut, bool FSR_in_Medium);
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
                int POI, double yabs_cut, bool FSR_in_Medium){
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
				c_entry.x0 = fourvec{0,x,y,0}; // initial position
				c_entry.x = c_entry.x0; // initial position
                c_entry.weight = weight;
                c_entry.is_vac = false;
				c_entry.is_virtual = false;
                c_entry.Tf = 0.0;
                c_entry.vcell.resize(3);
                c_entry.vcell[0] = 0.; 
                c_entry.vcell[1] = 0.; 
                c_entry.vcell[2] = 0.; 
				if (FSR_in_Medium){
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

		                    particle vp;
							vp.pid = 21;
							vp.mass = 0.0;
							vp.weight = weight;
							vp.p0 = k; 
							vp.p = vp.p0;
							vp.x0 = c_entry.x;
							vp.x = vp.x0;
							vp.T0 = 0.0; // not used
							vp.mfp0 = 0.0; // not used
							vp.is_vac = true;
							vp.is_virtual = true;

		                    c_entry.radlist.push_back(vp);
		                }
		            }
				}
                c_entry.p0 = p0; 
				c_entry.p = c_entry.p0; 
                plist.push_back(c_entry); 
            }
        }
      }

    // get the correct normalization for cross-section
    double normW = 0.5*pythia.info.sigmaGen()/pythia.info.weightSum();
    for (auto & p : plist) p.weight *= normW;
    LOG_INFO << "norm = " << normW;
}


#endif
