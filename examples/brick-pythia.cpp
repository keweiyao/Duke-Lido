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
#include "pythia_wrapper.h"

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
        std::vector<particle> plist;

        /// HardGen
        HardGen hardgen(
                args["pythia-setting"].as<fs::path>().string(), 
                args["ic"].as<fs::path>().string(),
                args["eid"].as<int>()           
        );
        hardgen.Generate(
                plist, 
                args["pythia-events"].as<int>(),
                4,
                2.5,
				true
        );
    
        /// Read Hydro
        Medium<2> med1(args["hydro"].as<fs::path>().string());

        /// Lido init
        initialize(table_mode,
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

        /// Assign each quark a transverse position according to TRENTo Nbin output
        /// Freestream particle to the start of hydro
        for (auto & p : plist){
            double vz2 = p.p.z()*p.p.z()/p.p.t()/p.p.t();
            double dt_fs = med1.get_tauH()/std::sqrt(1. - vz2);
            p.freestream(dt_fs);
        }



        // initial pT
        double pTi = mean_pT(plist);

        // run
        int counter = 0;
        int Ns = 10;

        while(med1.load_next()){
            double current_hydro_clock = med1.get_tauL();
            double hydro_dtau = med1.get_hydro_time_step();
            for (int i=0; i<Ns; ++i){
                double dtau = hydro_dtau/Ns; // use smaller dt step
                for (auto & p : plist){
                    if (p.freezeout) continue; // do not touch freezeout ones

                    // determine dt needed to evolve to the next tau
                    double tau = std::sqrt(p.x.t()*p.x.t()-p.x.z()*p.x.z());
                    double dt_lab = calcualte_dt_from_dtau(p.x, p.p, tau, dtau);

                    // get hydro information
                    double T = 0.0, vx = 0.0, vy = 0.0, vz = 0.;
                    med1.interpolate(p.x, T, vx, vy, vz);
                    double vabs = std::sqrt(vx*vx + vy*vy + vz*vy);
                    // regulate v
                    if (vabs > 1.){
                        LOG_WARNING << "regulate |v| = " << vabs << " > 1."; 
                        if (vabs > 1.-1e-6) {
                            double rescale = (1.-1e-6)/vabs;
                            vx *= rescale;
                            vy *= rescale;    
                            vz *= rescale;    
                        }
                    }
                    // x,p-update
                    int channel = 
                        OneBodyUpdate_Parton(dt_lab, T, {vx, vy, vz}, p);
                }
            }
            counter ++;
        }

        // final pT
        double pTf = mean_pT(plist);
        LOG_INFO << "Nparticles: " << plist.size();
        LOG_INFO << "Initial pT: " << pTi << " GeV";
        LOG_INFO << "Final pT: " << pTf << " GeV";

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



