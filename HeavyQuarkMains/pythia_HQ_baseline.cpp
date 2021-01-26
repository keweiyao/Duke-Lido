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
#include "lido.h"
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
            po::value<int>()->value_name("INT")->default_value(100,"100"),
           "number of Pythia events")
          ("ic,i",
            po::value<fs::path>()->value_name("PATH")->required(),
           "trento initial condition file")
          ("eid,j",
           po::value<int>()->value_name("INT")->default_value(0,"0"),
           "trento event id")
           ("output,o",
           po::value<fs::path>()->value_name("PATH")->default_value("./"),
           "output file prefix or folder")
	   ("Q0,q",
           po::value<double>()->value_name("DOUBLE")->default_value(.5,".5"),
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
            if (!fs::exists(args["pythia-setting"].as<fs::path>())) {
                throw po::error{"<pythia-setting> path does not exist"};
                return 1;
            }
        }

        std::vector<double> TriggerBin({
           2,3,4,6,8,10,12,14,16,18,20,22,24,26,28,30,
           35,40,45,50,55,60,65,70,75,80,85,90,100,
           110,120,140,160,180,200,
           240,280,320,360,400,500});

        double Q0 = args["Q0"].as<double>();
                
        // Fill in all events
        std::vector<particle> plist, dlist;
        for (int iBin = 0; iBin < TriggerBin.size()-1; iBin++) {
            HQGenerator pythiagen(
                            args["pythia-setting"].as<fs::path>().string(),
                            args["ic"].as<fs::path>().string(),
                            TriggerBin[iBin],
                            TriggerBin[iBin+1],
                            args["eid"].as<int>(),
                            Q0
                            );
            dlist.clear();
            pythiagen.Generate(dlist,
                             args["pythia-events"].as<int>(),
                             4.0);
            plist.insert(plist.end(), dlist.begin(), dlist.end());
        }

	for (auto & p : plist) {
            double tau = p.x.x0();
            double etas = p.x.x3();
            p.x.a[0] = tau*std::cosh(etas);
            p.x.a[3] = tau*std::sinh(etas);
            double vzcell = std::tanh(etas);
            double gamma = 1./std::sqrt(1.-vzcell*vzcell);
            double vz0cell = std::tanh(p.x0.x3());
            p.p = p.p.boost_back(0,0,vzcell);
            p.p0 = p.p0.boost_back(0,0,vz0cell);
            p.vcell[0] = p.vcell[0]/gamma/(1+vzcell*p.vcell[2]);
            p.vcell[1] = p.vcell[1]/gamma/(1+vzcell*p.vcell[2]);     
            p.vcell[2] = (vzcell+p.vcell[2])/(1+vzcell*p.vcell[2]);  

        }

        int processid = getpid();
        std::stringstream outputfilename1,outputfilename2;
        outputfilename1 << args["output"].as<fs::path>().string() 
                        << "c-quark-frzout.dat";
        outputfilename2 << args["output"].as<fs::path>().string()
                        << "b-quark-frzout.dat";
        output_oscar(plist, 4, outputfilename1.str());
        output_oscar(plist, 5, outputfilename2.str());
  
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



