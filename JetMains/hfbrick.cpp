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
#include "pythia_jet_gen.h"
#include "Hadronize.h"
#include "jet_finding.h"

namespace po = boost::program_options;
namespace fs = boost::filesystem;

void output_jet(std::string fname, std::vector<particle> plist){
    std::ofstream f(fname);
    for (auto & p : plist) f << p.pid << " " << p.p << " " << p.x << " " << p.p0 << std::endl;
    f.close();
}

struct event{
std::vector<particle> plist, colorlist, thermal_list, hlist;
std::vector<current> clist;
std::vector<HadronizeCurrent> slist;
double sigma;
};

int main(int argc, char* argv[]){
    using OptDesc = po::options_description;
    OptDesc options{};
    options.add_options()
          ("help", "show this help message and exit")
          ("nevents,n",
            po::value<int>()->value_name("INT")->default_value(100,"100"),
           "number of Pythia events")
          ("lido-setting,s", 
            po::value<fs::path>()->value_name("PATH")->required(),
           "Lido table setting file")
          ("lido-table,t", 
            po::value<fs::path>()->value_name("PATH")->required(),
           "Lido table path to file")  
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
     
        /// use process id to define filename
        int processid = getpid();
        std::stringstream fheader;
        fheader << args["output"].as<fs::path>().string() 
                << processid << "";
        /// Initialize Lido in-medium transport
        initialize(table_mode,
            args["lido-setting"].as<fs::path>().string(),
            args["lido-table"].as<fs::path>().string()
        );
        
        
        std::vector<particle> plist, pOut_list;
        fourvec x0{0,0,0,0}, p0{30,0,0,std::sqrt(30*30-1.3*1.3)};

        for (int i=0; i<args["nevents"].as<int>(); i++){
            particle p;
            p.p = p0;
            p.p0 = p0;
            p.x = x0;
            p.x0 = x0;
            p.tau_i = 0;
            p.pid = 4;
            p.mass = 1.3;
            p.col = 101;
            p.acol = 0;
            p.vcell.resize(3);
            p.is_virtual=false;
            plist.push_back(p);
        }
        double dt = 0.05*5.076;
        for (int i=0; i<101; i++){
            //LOG_INFO << "S";
            double t = i*dt;
            LOG_INFO << i;
            for (auto & p : plist){  
                //LOG_INFO << p.p;   
                // get hydro information
                double T = .25, vx = 0.0, vy = 0.0, vz = 0.0;
                int fs_size = update_particle_momentum_Lido(
                            dt, T, {vx, vy, vz}, p, pOut_list);     
            }
            
            if (i%20==0){
                LOG_INFO << "h";
                std::stringstream fname;
                fname << fheader.str()<<"-"<< (i/20) <<".dat";
                output_jet(fname.str(), plist);
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



