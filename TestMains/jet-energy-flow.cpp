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
#include <unistd.h>

#include "simpleLogger.h"
#include "Medium_Reader.h"
#include "workflow.h"
#include "PGunWithShower.h"
#include "Hadronize.h"
#include "jet_finding.h"

namespace po = boost::program_options;
namespace fs = boost::filesystem;

void output_jet(std::string fname, std::vector<particle> plist){
    std::ofstream f(fname);
    for (auto & p : plist) f << p.pid << " " << p.p << " " << p.x << " " << p.Q0 << std::endl;
    f.close();
}

int main(int argc, char* argv[]){
    using OptDesc = po::options_description;
    OptDesc options{};
    options.add_options()
          ("help", "show this help message and exit")
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
           ("output,o",
           po::value<fs::path>()->value_name("PATH")->default_value("./"),
           "output file prefix or folder")
          ("Q0,q",
           po::value<double>()->value_name("DOUBLE")->default_value(.4,".4"),
           "Scale [GeV] to insert in-medium transport")
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

        /// use process id to define filename
        int processid = getpid();
        /// Initialize Lido in-medium transport
        initialize(table_mode,
            args["lido-setting"].as<fs::path>().string(),
            args["lido-table"].as<fs::path>().string()
        );

        // Scale to insert In medium transport
        std::cout << "set Q0" << std::endl;
        double Q0 = args["Q0"].as<double>();
        PGunWShower HardGen(Q0);

    for (int ie=0; ie<args["pythia-events"].as<int>(); ie++){
        std::vector<particle> plist, hlist, thermal_list,
                              new_plist, pOut_list;
        std::vector<current> clist;
        std::vector<HadronizeCurrent> slist;
        // Initialize parton list from pythia
        HardGen.Generate(100, plist);
        /// Initialzie a hydro reader
        Medium<2> med1(args["hydro"].as<fs::path>().string());
        // freestream form t=0 to tau=tau0
        for (auto & p : plist){
            p.Tf = 0.16;
            fourvec x0{0.,0.,0.,0.};
            p.x0 = x0;
            p.x = x0;
            if (p.x.tau() < med1.get_tauH()){
                p.freestream(
                    compute_realtime_to_propagate(
                        med1.get_tauH(), p.x, p.p
                    )
                );
            }
        }   
        int Nstep = 0, iFrame=0;
        while(med1.load_next()){
            double current_hydro_clock = med1.get_tauL();
            double hydro_dtau = med1.get_hydro_time_step();
            double dtau = hydro_dtau, DeltaTau;
            new_plist.clear();
            for (auto & p : plist){     
                if (p.Tf < 0.15) {
                    new_plist.push_back(p);
                    continue;       
                }
                if (p.x.tau() > current_hydro_clock+dtau){
                    new_plist.push_back(p);
                    continue;
                }
                else{
                    DeltaTau = current_hydro_clock 
                             + dtau - p.x.tau();
                }
                // get hydro information
                double T = 0.0, vx = 0.0, vy = 0.0, vz = 0.0;
                med1.interpolate(p.x, T, vx, vy, vz);
                double vzgrid = p.x.z()/p.x.t();
                fourvec ploss = p.p;
                int fs_size = update_particle_momentum_Lido(
                      DeltaTau, T, {vx, vy, vz}, p, pOut_list);
                 
                for (auto & fp : pOut_list) {
                    ploss = ploss - fp.p;
                    new_plist.push_back(fp);
                }   
            }
            if (Nstep%10==0){
                std::stringstream fname;
                fname << args["output"].as<fs::path>().string()
                          << "/" << iFrame << "/"
                          << processid << "-" << ie << ".dat";
                output_jet(fname.str(), plist);
                iFrame ++;
            }
            plist = new_plist;
            Nstep ++;
         }
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



