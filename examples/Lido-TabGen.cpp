#include <string>
#include <iostream>
#include <exception>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/program_options.hpp>

#include "simpleLogger.h"
#include "workflow.h"

namespace po = boost::program_options;
namespace fs = boost::filesystem;

int main(int argc, char* argv[]){
    using OptDesc = po::options_description;
    OptDesc options{};
    options.add_options()
          ("help,h", "show this help message and exit")
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
        if (!args.count("lido-setting")){
            throw po::required_option{"<lido-setting>"};
            return 1;
        }
        else {
            if (!fs::exists(args["lido-setting"].as<fs::path>())){
                throw po::error{"<lido-setting> path does not exist"};
                return 1;
            }
        }
        if (!args.count("lido-table")){
            throw po::required_option{"<lido-table>"};
            return 1;
        }

        initialize("new", 
            args["lido-setting"].as<fs::path>().string(), 
            args["lido-table"].as<fs::path>().string(), 
            args["mu"].as<double>(),
            args["afix"].as<double>(), 
            args["k-factor"].as<double>(), 
            args["t-scale"].as<double>(), 
            args["e-scale"].as<double>(), 
            args["t-power"].as<double>(), 
            args["e-power"].as<double>(), 
            args["gamma"].as<double>(), 
            args["qcut"].as<double>(), 
            args["rvac"].as<double>()
        );
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


