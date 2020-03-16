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


#include "simpleLogger.h"
#include "Medium_Reader.h"
#include "workflow.h"
#include "pythia_HQ_gen.h"

namespace po = boost::program_options;
namespace fs = boost::filesystem;

int main(int argc, char* argv[]){
    using OptDesc = po::options_description;
    OptDesc options{};
    options.add_options()
          ("help", "show this help message and exit")
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
          ("hydro",
           po::value<fs::path>()->value_name("PATH")->required(),
           "hydro file")
          ("lido-setting,s", 
            po::value<fs::path>()->value_name("PATH")->required(),
           "Lido table setting file")
          ("lido-table,t", 
            po::value<fs::path>()->value_name("PATH")->required(),
           "Lido table path to file")   
	  ("cl",
           po::value<int>()->value_name("INT")->default_value(0,"0"),
           "cent low")
	  ("ch",
           po::value<int>()->value_name("INT")->default_value(0,"100"),
           "cent high")
           ("output,o",
           po::value<fs::path>()->value_name("PATH")->default_value("./"),
           "output file prefix or folder")
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
            table_mode = (fs::exists(args["lido-table"].as<fs::path>())) ? 
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
        std::vector<particle> plist, dlist, flist;

    
        /// Read Hydro
        Medium<2> med1(args["hydro"].as<fs::path>().string());

        /// Lido init
        initialize(table_mode,
                args["lido-setting"].as<fs::path>().string(),
                args["lido-table"].as<fs::path>().string()
                );

        /// generate events from pythia
        std::vector<double> pThatbins({2,4,6,8,10,15,20,30,40,50,60,80,100,120,150,200,300,500});
        for (int i=0; i<pThatbins.size()-1; i++){
            HQGenerator hardgen(args["pythia-setting"].as<fs::path>().string(),
                            args["ic"].as<fs::path>().string(),
                            args["eid"].as<int>(),
                            pThatbins[i], pThatbins[i+1]
                            );
            dlist.clear();
            hardgen.Generate(dlist,
                             args["pythia-events"].as<int>(),
                             3.0);
            plist.insert(plist.end(), dlist.begin(), dlist.end());
        }

        // freestream form t=0 to tau=tau0
        for (auto & p : plist){
            p.Tf = 0.161;
            p.origin=-1;
            if (p.x.tau() < med1.get_tauH()){
                p.freestream(
                    compute_realtime_to_propagate(
                        med1.get_tauH(), p.x, p.p
                    )
                );
            }
        }   
        int itau = 0;
        while(med1.load_next()){
            
            double current_hydro_clock = med1.get_tauL();
            double hydro_dtau = med1.get_hydro_time_step();
            //LOG_INFO << current_hydro_clock/fmc_to_GeV_m1 
            //        << " [fm/c]\t" 
            //       << " # of hard =" << plist.size();
            // further divide hydro step into 10 transpor steps
            int Ns = 2; 
            double dtau = hydro_dtau/Ns, DeltaTau;
            for (int i=0; i<Ns; ++i){
                for (auto & p : plist){     
                    if (p.Tf < 0.161) {
                        continue;       
                    }
                    if (p.x.tau() > current_hydro_clock+(i+1)*dtau){
                        // if the particle time is in the future 
                        // (not formed to the medium yet)
                        continue;
                    }
                    else{
                        DeltaTau = current_hydro_clock 
                                   + (i+1)*dtau - p.x.tau();
                    }
                    // get hydro information
                    double T = 0.0, vx = 0.0, vy = 0.0, vz = 0.0;
                    double vzgrid = p.x.z()/p.x.t();
                    med1.interpolate(p.x, T, vx, vy, vz);

                    int fs_size = update_particle_momentum_Lido(
                                  DeltaTau, T, {vx, vy, vz}, p, flist);
                }
            }
            itau ++; 
        }
        int processid = getpid();
        std::stringstream outputfilename1;
        outputfilename1 << args["output"].as<fs::path>().string()
                        << "c-quark-" << processid <<".dat";
        output_oscar(plist, 4, outputfilename1.str());
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



