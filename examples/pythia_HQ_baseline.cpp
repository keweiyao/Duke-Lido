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
#include <unistd.h>

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
        std::vector<particle> plist, dlist;
        plist.clear();
        /// HardGen
        std::vector<double> pThatbins({2,4,6,8,10,15,20,30,40,50,60,80,100,120,150,200,300});
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
        int processid = getpid();
        std::stringstream outputfilename1,outputfilename2;
        outputfilename1 << "c-quark-" << processid<<".dat";
        outputfilename2 << "b-quark-" << processid<<".dat";
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



