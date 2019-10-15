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
#include <sstream>

#include "simpleLogger.h"
#include "Medium_Reader.h"
#include "workflow.h"
#include "pythia_jet_gen.h"

namespace po = boost::program_options;
namespace fs = boost::filesystem;

void output_jet(std::string fname, std::vector<particle> plist){
    std::ofstream f(fname);
    for (auto & p : plist) f << p.pid << " " << p.p << " " << p.weight << std::endl;
    f.close();
}

int main(int argc, char* argv[]){
    using OptDesc = po::options_description;
    OptDesc options{};
    options.add_options()
          ("help", "show this help message and exit")
          ("pthat-low,l",
            po::value<int>()->value_name("INT")->default_value(100,"100"),
           "pThatMin")
          ("pthat-high,h",
            po::value<int>()->value_name("INT")->default_value(120,"120"),
           "pThatMax")
          ("pythia-setting,y",
           po::value<fs::path>()->value_name("PATH")->required(),
           "Pythia setting file")
          ("pythia-events,n",
            po::value<int>()->value_name("INT")->default_value(100,"100"),
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

        /// HardGen
        LOG_INFO<<"hard gen";
        PythiaGen pythiagen(
                args["pythia-setting"].as<fs::path>().string(),
                args["ic"].as<fs::path>().string(),
                args["pthat-low"].as<int>(),
                args["pthat-high"].as<int>(),
                args["eid"].as<int>()
        );

        /// Lido init
        LOG_INFO<<"init lido";
        initialize(table_mode,
                args["lido-setting"].as<fs::path>().string(),
                args["lido-table"].as<fs::path>().string()
                );

        int Ns = 10;
        for (int ie=0; ie<args["pythia-events"].as<int>(); ie++){
            std::vector<particle> plist, new_plist, pOut_list;
            pythiagen.Generate(plist);
            /// Read Hydro
            Medium<2> med1(args["hydro"].as<fs::path>().string());
            // freestream form t=0 to tau=tau0
            for (auto & p : plist) {
              p.freestream(compute_realtime_to_propagate(med1.get_tauH(), p.x, p.p));
              p.Tf = 0.151;
            }
            while(med1.load_next()){
                double current_hydro_clock = med1.get_tauL();
                double hydro_dtau = med1.get_hydro_time_step();
                double dtau = hydro_dtau/Ns;
               // LOG_INFO << current_hydro_clock/5.026 << " [fm/c]\t" 
                //         << " #=" << plist.size();
                for (int i=0; i<Ns; ++i){
                    new_plist.clear();
                    for (auto & p : plist){
                        // get hydro information
                        double T = 0.0, vx = 0.0, vy = 0.0, vz = 0.0;
                        med1.interpolate(p.x, T, vx, vy, vz);
                        double vabs = std::sqrt(vx*vx + vy*vy + vz*vy);
                        // regulate v
                        if (vabs > 1.-1e-6){
                           // LOG_WARNING << "regulate |v| = " 
                             //           << vabs << " > 1. and "
                              //          << "y = " 
                        //<< 0.5*std::log((p.x.t()+p.x.z())/(p.x.t()-p.x.z()));
                            double rescale = (1.-1e-6)/vabs;
                            vx *= rescale;
                            vy *= rescale;    
                            vz *= rescale;                            
                        }
                        int fs_size = update_particle_momentum_Lido(
                                  dtau, T, {vx, vy, vz}, p, pOut_list);
                        for (auto & k : pOut_list) new_plist.push_back(k);
                    }
                    plist = new_plist;
                }
            }
            std::ostringstream s;
            s << "results/" << args["pthat-low"].as<int>() << "-"
              << args["pthat-high"].as<int>() << "-" << ie;
            output_jet(s.str(), plist);
        }
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



