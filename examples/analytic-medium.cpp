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
#include "workflow.h"

namespace po = boost::program_options;
namespace fs = boost::filesystem;

void output_jet(std::string fname, std::vector<particle> plist){
    std::ofstream f;
    f.open(fname, std::ofstream::out | std::ofstream::app);
    f << "# " << plist.size() << std::endl;
    for (auto & p : plist) f << p.pid << " " << p.x << " " << p.p << " " << p.origin << std::endl;
    f.close();
}

int main(int argc, char* argv[]){
    using OptDesc = po::options_description;
    OptDesc options{};
    options.add_options()
          ("help", "show this help message and exit")
          ("nparticles",
            po::value<int>()->value_name("INT")->default_value(100,"100"),
           "number of hard partons")
          ("nevents",
            po::value<int>()->value_name("INT")->default_value(10,"10"),
           "number of events")
          ("energy",
            po::value<double>()->value_name("FLOAT")->default_value(10.0,"10.0"),
           "fixed parton energy [GeV]")
          ("temp",
           po::value<double>()->value_name("FLOAT")->default_value(0.3,"0.3"),
           "medium temperature at t=t0 [GeV]")
          ("taui",
           po::value<double>()->value_name("FLOAT")->default_value(0.0,"0.0"),
           "medium starting time [fm/c]")
          ("tauf",
           po::value<double>()->value_name("FLOAT")->default_value(5.0,"5.0"),
           "medium stopping time [fm/c]")
          ("dtau",
           po::value<double>()->value_name("FLOAT")->default_value(0.1,"0.1"),
           "hydro time step dtau [fm/c]")
          ("nu",
           po::value<double>()->value_name("FLOAT")->default_value(0.5,"0.5"),
           "T = T0 * (tau0/tau) ^ (2-1/nu)")
          ("lido-setting,s",
            po::value<fs::path>()->value_name("PATH")->required(),
           "Lido table setting file")
          ("lido-table,t",
            po::value<fs::path>()->value_name("PATH")->required(),
           "Lido table path to file")
          ("mass",
            po::value<double>()->value_name("FLOAT")->default_value(1.3,"1.3"),
            "quark mass")
          ("pid",
            po::value<int>()->value_name("INT")->default_value(4,"4"),
            "particle pid, 4 or 21")
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
        // check pid
        if (!args.count("pid")){
            if ( (std::abs(args["pid"].as<int>()) != 4)
                && (std::abs(args["pid"].as<int>()) != 21) ){
                throw po::error{"pid only suppots 4, 21"};
                return 1;
            }
        }

        /// Lido init
        initialize(table_mode,
                args["lido-setting"].as<fs::path>().string(),
                args["lido-table"].as<fs::path>().string()
                );

        /// Assign each quark a transverse position according to TRENTo Nbin output
        /// Freestream particle to the start of hydro
        double ti = args["taui"].as<double>() * 5.026; // convert to GeV^-1
        double tf = args["tauf"].as<double>() * 5.026; // convert to GeV^-1
        double dt = args["dtau"].as<double>() * 5.026; // convert to GeV^-1
        double T0 = args["temp"].as<double>();
        double nu = args["nu"].as<double>();
        int Nsteps = int((tf-ti)/dt)+1;


        std::vector<particle> plist, new_plist, temp_list;
        int pid = args["pid"].as<int>();
        double E0 = args["energy"].as<double>();
        double M = args["mass"].as<double>();
        fourvec p0{E0, 0, 0, std::sqrt(E0*E0-M*M)},
                x0{0,0,0,0};

        std::ofstream sourcef;
        sourcef.open("results/single_recoil.dat",
                      std::ofstream::out | std::ofstream::app);
        // start
        
        for (int ie=0; ie<args["nevents"].as<int>(); ie++){
            plist.clear();
            double deltaE = 0;
            for (int i=0; i<args["nparticles"].as<int>(); i++){
                particle p_init;
                p_init.origin = -1;
                p_init.p = p0; //
                p_init.p0 = p0; //
                p_init.mass = M; //
                p_init.pid = pid;
                p_init.x0 = x0; //
                p_init.x = x0; //
                p_init.weight = 1.;
                p_init.Tf = 0.0;
                p_init.is_vac = false;
                p_init.is_virtual = false;
                p_init.vcell.resize(3);
                p_init.vcell[0] = 0.;
                p_init.vcell[1] = 0.;
                p_init.vcell[2] = 0.;
                plist.push_back(p_init);
            }

            for (auto & p : plist) p.freestream(compute_realtime_to_propagate(ti, p.x, p.p));

            for (int i=0; i<Nsteps; i++){
                new_plist.clear();
                double t = ti + dt*i;
                double T = T0 * std::pow(ti/t, 2./3.-1./nu/3.);
                if (i%100==0) LOG_INFO << t/5.026 << " [fm/c]\t"<< T
                                       << " [GeV]" << " #=" << plist.size();
                for (auto & p : plist){
                    p.p = p.p0;
                    int n = update_particle_momentum_Lido(
                                  dt, T, {0., 0., 0.}, p, temp_list);
                    for (auto & k : temp_list) {
                        if (k.origin == -1) {
                            if(t>=(ti+tf)/2.)
                                deltaE += 
                                     (E0 - k.p.t())
                                     /args["nparticles"].as<int>();
                            new_plist.push_back(k);
                        }
                    }
                }
                plist = new_plist;
            }

            std::ostringstream s;
            s << "results/test.dat";
            std::ofstream f;
            f.open(s.str(),  std::ofstream::out | std::ofstream::app);
            f << E0 << " " << (tf-ti)/2. << " " << deltaE << std::endl;
            f.close();
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
