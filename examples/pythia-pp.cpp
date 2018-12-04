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
                c_entry.freezeout = true;
                c_entry.Tf = 0.0;
                c_entry.vcell.resize(3);
                c_entry.vcell[0] = 0.; 
                c_entry.vcell[1] = 0.; 
                c_entry.vcell[2] = 0.; 
                c_entry.p = p0; // without FSR energy loss
                
                plist.push_back(c_entry); 
            }
        }
      }

    // get the correct normalization for cross-section
    double normW = 0.5*pythia.info.sigmaGen()/pythia.info.weightSum();
    for (auto & p : plist) p.weight *= normW;
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
    ;

    po::variables_map args{};
    try{
        po::store(po::command_line_parser(argc, argv).options(options).run(), args);
        if (args.count("help")){
                std::cout << "usage: " << argv[0] << " [options]\n"
                          << options;
                return 0;
        }    
        // check trento event
        if (!args.count("ic")){
            throw po::required_option{"<ic>"};
            return 1;
        }
        else{
            if (!fs::exists(args["ic"].as<fs::path>())){
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
            if (!fs::exists(args["pythia-setting"].as<fs::path>())){
                throw po::error{"<pythia-setting> path does not exist"};
                return 1;
            }
        }

        // start
        std::vector<particle> plist;

        /// HardGen
        HardGen hardgen(args["pythia-setting"].as<fs::path>().string(),
                            args["ic"].as<fs::path>().string(),
                            args["eid"].as<int>()
                            );
        hardgen.Generate(plist,
                            args["pythia-events"].as<int>(),
                            4,
                            2.5);


        
        double pTf = mean_pT(plist);
        LOG_INFO << "Nparticles: " << plist.size();
        LOG_INFO << "Mean pT: " << pTf << " GeV";

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



