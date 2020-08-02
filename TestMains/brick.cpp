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
#include "workflow.h"

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
          ("lido-setting,s", 
            po::value<fs::path>()->value_name("PATH")->required(),
           "Lido table setting file")
          ("lido-table,t", 
            po::value<fs::path>()->value_name("PATH")->required(),
           "Lido table path to file")  
           ("output,o",
           po::value<fs::path>()->value_name("PATH")->default_value("./"),
           "output file prefix or folder")
           ("einit,e",
           po::value<double>()->value_name("DOUBLE")->default_value(100.),
           "initial parton energy")
           ("temp",
           po::value<double>()->value_name("DOUBLE")->default_value(.3),
           "medium temperature")
           ("pid,i",
           po::value<int>()->value_name("INT")->default_value(21),
           "initial parton pid")
	  ("muT",
           po::value<double>()->value_name("DOUBLE")->default_value(1.5,"1.5"),
	   "mu_min/piT")
	  ("afix",
           po::value<double>()->value_name("DOUBLE")->default_value(-1.,"-1."),
           "fixed alpha_s, <0 for running alphas")
           ("cut",
           po::value<double>()->value_name("DOUBLE")->default_value(4.,"4."),
           "cut between diffusion and scattering, Qc^2 = cut*mD^2")
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

        /// use process id to define filename
        int processid = getpid();
        /// Initialize Lido in-medium transport
        double muT = args["muT"].as<double>();
        double cut = args["cut"].as<double>();
        double afix = args["afix"].as<double>();

        initialize(table_mode, 
            args["lido-setting"].as<fs::path>().string(), 
            args["lido-table"].as<fs::path>().string(), 
	    muT, afix, 4, cut);

        int pid = args["pid"].as<int>();
        double T0 = args["temp"].as<double>();
        double E0 = args["einit"].as<double>();
        double L0 = 100*5.076;
        double dL = L0/5000;
        std::vector<particle> plist, new_plist, pOut_list;
        plist.resize(1000);
        fourvec x0{0.,0.,0.,0.};
        fourvec p0{E0,0.,0.,E0};

        for (auto & p : plist){
            p.pid = pid;
            p.T0 = T0;
            p.Tf = T0;
            p.x0 = x0;
            p.x = x0;
            p.tau_i = 0.;
            p.p0 = p0;
            p.p = p0;
            p.vcell.resize(3);
            p.is_virtual = false;
        }
        for (double l=0; l<L0; l+=dL){
            LOG_INFO << l;
            for (auto & p : plist){
                p.p = p0;
                int fs_size = update_particle_momentum_Lido(
                      dL, T0, {0,0,0}, p, pOut_list);       
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



