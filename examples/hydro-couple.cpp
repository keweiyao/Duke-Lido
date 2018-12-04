#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <exception>
#include <random>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/program_options.hpp>

#include "Pythia8/Pythia.h"
#include "simpleLogger.h"
#include "Medium_Reader.h"
#include "workflow.h"

namespace po = boost::program_options;
namespace fs = boost::filesystem;

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
                c_entry.p0 = p0; // 
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
                
                plist.push_back(c_entry); 
            }
        }
      }

    // get the correct normalization for cross-section
    double normW = 0.5*pythia.info.sigmaGen()/pythia.info.weightSum();
    for (auto & p : plist) p.weight *= normW;
    LOG_INFO << "norm = " << normW;
}

bool is_file_exist(std::string fileName)
{
    std::ifstream infile(fileName);
    return infile.good();
}

int main(int argc, char* argv[]){
    using OptDesc = po::options_description;
    OptDesc options{};
    options.add_options()
          ("help,h", "show this help message and exit")
          ("pythia-setting,y",
           po::value<fs::path>()->value_name("PATH")->required(),
           "Pythia setting file")
          ("pythia-events,n",
            po::value<int>()->value_name("INT")->default_value(50000,"50000"),
           "number of Pythia events")
          ("ic,i",
            po::value<fs::path>()->value_name("PATH")->required(),
           "trento initial condition file")
          ("eid,j",
           po::value<int>()->value_name("INT")->default_value(0,"0"),
           "trento event id")
          ("hydro,h",
           po::value<fs::path>()->value_name("PATH")->required(),
           "hydro file")
          ("lido-setting,s", 
            po::value<fs::path>()->value_name("PATH")->required(),
           "Lido table setting file")
          ("lido-table,t", 
            po::value<fs::path>()->value_name("PATH")->required(),
           "Lido table path to file")
          ("mu,m", 
            po::value<double>()->value_name("FLOAT")->default_value(1.0,"1.0"),
            "medium scale paramtmer")
          ("afix,f", 
            po::value<double>()->value_name("FLOAT")->default_value(-1.0,"-1.0"),
            "fixed alphas value, -1 is running")
          ("k-factor,k", 
            po::value<double>()->value_name("FLOAT")->default_value(0.0,"0.0"),
            "K-factor for the delta-qhat")
          ("t-scale,a", 
            po::value<double>()->value_name("FLOAT")->default_value(1.0,"1.0"),
            "rescale the T-dependence")
          ("e-scale,b", 
            po::value<double>()->value_name("FLOAT")->default_value(1.0,"1.0"),
            "rescale the p-dependence")
          ("t-power,p", 
            po::value<double>()->value_name("FLOAT")->default_value(1.0,"1.0"),
            "T-dependence power")
          ("e-power,q", 
            po::value<double>()->value_name("FLOAT")->default_value(1.0,"1.0"),
            "p-dependence power")
          ("gamma,g", 
            po::value<double>()->value_name("FLOAT")->default_value(0.0,"0.0"),
            "kpara / kperp anisotropy parameter")
          ("qcut,c",
            po::value<double>()->value_name("FLOAT")->default_value(1.0,"1.0"), 
            "separation scale Q2 = qcut * mD2")
          ("rvac,r", 
            po::value<double>()->value_name("FLOAT")->default_value(1.0,"1.0"), 
            "vacuum-like radiation remove factor")    
    ;
    po::variables_map args{};
    try{
        po::store(po::command_line_parser(argc, argv).options(options).run(), args);
        if (args.count("help")){
                std::cout << "usage: " << argv[0] << " [options]\n"
                          << options;
                return 0;
        }    
        // check lido setting
        if (!args.count("lido-setting")){
            throw po::required_option{"<lido-setting>"};
            return 1;
        }
        else{
            if (!fs::exists(args["lido-setting"].as<fs::path>())){
                throw po::error{"<lido-setting> path does not exist"};
                return 1;
            }
        }
        // check lido table
        std::string table_mode;
        if (!args.count("lido-table")){
            throw po::required_option{"<lido-table>"};
            return 1;
        }
        else{
            table_mode = (fs::exists(args["lido-setting"].as<fs::path>())) ? 
                         "old" : "new";
        }
        // check trento event
        if (!args.count("ic")){
            throw po::required_option{"<ic>"};
            return 1;
        }
        else{
            if (!fs::exists(args["ic"].as<fs::path>())) {
                throw po::error{"<ic> path does not exist"};
                return 1;
            }
        }
        // check pythia setting
        if (!args.count("pythia-setting")){
            throw po::required_option{"<pythia-setting>"};
            return 1;
        }
        else{
            if (!fs::exists(args["pythia-setting"].as<fs::path>())) {
                throw po::error{"<pythia-setting> path does not exist"};
                return 1;
            }
        }
        // check hydro setting
        if (!args.count("hydro")){
            throw po::required_option{"<hydro>"};
            return 1;
        }
        else{
            if (!fs::exists(args["hydro"].as<fs::path>())) {
                throw po::error{"<hydro> path does not exist"};
                return 1;
            }
        }

        // start
        std::vector<particle> plist;

        /// HardGen
        HardGen hardgen(
                args["pythia-setting"].as<fs::path>().string(), 
                args["ic"].as<fs::path>().string(),
                args["eid"].as<int>()           
        );
        hardgen.Generate(
                plist, 
                args["pythia-events"].as<int>(),
                4,
                2.5
        );
    
        /// Read Hydro
        Medium<2> med1(args["hydro"].as<fs::path>().string());
        
        /// Assign each quark a transverse position according to TRENTo Nbin output
        /// Freestream particle to the start of hydro
        for (auto & p : plist){
            double vz2 = p.p.z()*p.p.z()/p.p.t()/p.p.t();
            double dt_fs = med1.get_tauH()/std::sqrt(1. - vz2);
            p.freestream(dt_fs);
        }

        /// Lido init
        initialize(table_mode,
                args["lido-setting"].as<fs::path>().string(),
                args["lido-table"].as<fs::path>().string(),
                args["mu"].as<double>(),
                args["afix"].as<double>(),
                args["k-factor"].as<double>(),
                args["t-scale"].as<double>(),
                args["p-scale"].as<double>(),
                args["t-power"].as<double>(),
                args["p-power"].as<double>(),
                args["gamma"].as<double>(),
                args["qcut"].as<double>(),
                args["rvac"].as<double>()
                );

        // initial pT
        double pTi = mean_pT(plist);

        // run
        int counter = 0;
        int Ns = 10;

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
                    // x-update
                    p.freestream(dt_lab);

                    // get hydro information
                    double T = 0.0, vx = 0.0, vy = 0.0, vz = 0.;
                    med1.interpolate(p.x, T, vx, vy, vz);
                    double vabs = std::sqrt(vx*vx + vy*vy + vz*vy);
                    // regulate v
                    if (vabs > 1.){
                        LOG_WARNING << "regulate |v| = " << vabs << " > 1."; 
                        if (vabs > 1.-1e-6) {
                            double rescale = (1.-1e-6)/vabs;
                            vx *= rescale;
                            vy *= rescale;    
                            vz *= rescale;    
                        }
                    }
                    // p-update
                    int channel = 
                        update_particle_momentum_Lido(dt_lab, T, {vx, vy, vz}, p);
                }
            }
            counter ++;
        }

        // final pT
        double pTf = mean_pT(plist);
        LOG_INFO << "Nparticles: " << plist.size();
        LOG_INFO << "Initial pT: " << pTi << " GeV";
        LOG_INFO << "Final pT: " << pTf << " GeV";

        output_oscar(plist, "c-quark-frzout.dat");
    }
    catch (const po::required_option& e){
        std::cout << e.what() << "\n";
        std::cout << "usage: " << argv[0] << " [options]\n"
                  << options;
        return 1;
    }
    catch (const std::exception& e) {
       // For all other exceptions just output the error message.
       std::cerr << e.what() << '\n';
       return 1;
    }    
    return 0;
}



