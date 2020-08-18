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
#include "Hadronize.h"
#include "jet_finding.h"

namespace po = boost::program_options;
namespace fs = boost::filesystem;

void output_jet(std::string fname, std::vector<particle> plist){
    std::ofstream f(fname);
    for (auto & p : plist) f << p.pid << " " << p.p << " " << p.x0 << " " << p.x << std::endl;
    f.close();
}

int main(int argc, char* argv[]){
    using OptDesc = po::options_description;
    OptDesc options{};
    options.add_options()
          ("help", "show this help message and exit")         
          ("response-table,r", 
            po::value<fs::path>()->value_name("PATH")->required(),
           "response table path to file")  
    ;
    po::variables_map args{};
    try{
        po::store(po::command_line_parser(argc, argv).options(options).run(), args);
        if (args.count("help")){
                std::cout << "usage: " << argv[0] << " [options]\n"
                          << options;
                return 0;
        }    
       
        // check lido-response table
        bool need_response_table;
        if (!args.count("response-table")){
            throw po::required_option{"<response-table>"};
            return 1;
        }
        else{
            need_response_table = (fs::exists(args["response-table"].as<fs::path>())) ? 
                         false : true;
        }

        /// use process id to define filename
        int processid = getpid();
        std::stringstream outputfilename2;

        std::vector<current> clist;
        auto ReDistributer = MediumResponse("Gmu");
        if (need_response_table) ReDistributer.init(args["response-table"].as<fs::path>().string());
        else ReDistributer.load(args["response-table"].as<fs::path>().string());

        double xs=0, ys=0, etas=0; 
        fourvec total{0,0,0,0};
                double T = 0.0, vx = 0.0, vy = 0.0, vz = 0.0;
                current J; 
                fourvec ploss{1,1,-1,1};
                ploss = ploss;//*(hydro_dtau/0.4/5.076);
                total = total + ploss;
                //ploss = ploss.boost_to(0, 0, vzgrid);
                J.p = ploss;              
                J.chetas = 1.;
                J.shetas = 0.;
                J.cs = std::sqrt(.333);
                clist.push_back(J);  
        std::stringstream outputfilename1;
        outputfilename1 << "gT.dat";
        std::vector<HadronizeCurrent> slist;
        slist.clear();
        TestSource(ReDistributer, clist, slist, outputfilename1.str());
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



