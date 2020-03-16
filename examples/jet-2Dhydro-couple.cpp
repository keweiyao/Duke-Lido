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
    for (auto & p : plist) f << p.pid << " " << p.p << " " << p.x0 << " " << p.x << std::endl;
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
          ("heavy", 
            po::value<int>()->value_name("INT")->default_value(0,"0"),
           "require pythia event contains charm(4) or bottom(5) quark")
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

        if (args["heavy"].as<int>() >= 4){
            if (args["heavy"].as<int>() == 4)
                LOG_INFO << "Requires Pythia parton level has charm quark(s)";
            else if (args["heavy"].as<int>() == 5)
                LOG_INFO << "Requires Pythia parton level has bottom quark(s)";
            else {
                throw po::error{"Heavy trigger invalided"};
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
        /// Initialize a simple hadronizer
        JetDenseMediumHadronize Hadronizer;


        std::vector<double> TriggerBin({20,40,60,
90,110,130,140,160,180,200,250,300,350,400,500,700,1000,1500,2000});
        //std::vector<double> TriggerBin({30,60,100,200,300,
        //                                400,500,600,800,1500});
        //std::vector<double> TriggerBin({80,100,120,125,135,
        //                               145,160,180,200,250,300,1000});
        for (int iBin = 0; iBin < TriggerBin.size()-1; iBin++){
            /// Initialize a pythia generator for each pT trigger bin
            PythiaGen pythiagen(
                args["pythia-setting"].as<fs::path>().string(),
                args["ic"].as<fs::path>().string(),
                TriggerBin[iBin],
                TriggerBin[iBin+1],
                args["eid"].as<int>()
            );
            LOG_INFO << pythiagen.sigma_gen() << " " << TriggerBin[iBin] << " "<< TriggerBin[iBin];
            for (int ie=0; ie<args["pythia-events"].as<int>(); ie++){
                std::vector<particle> plist, hlist, thermal_list,
                                      new_plist, pOut_list;
                std::vector<current> clist;
                // Initialize parton list from python
                pythiagen.Generate(plist, args["heavy"].as<int>());
                double sigma_gen = pythiagen.sigma_gen()
                                  / args["pythia-events"].as<int>();
                

                /// Initialzie a hydro reader
                Medium<2> med1(args["hydro"].as<fs::path>().string());

                // freestream form t=0 to tau=tau0
                for (auto & p : plist){
                    p.Tf = 0.17;
                    p.origin=-1;
                    if (p.x.tau() < med1.get_tauH()){
                        p.freestream(
                            compute_realtime_to_propagate(
                                med1.get_tauH(), p.x, p.p
                            )
                        );
                    }
                }   

                // Energy-momentum checkbook
                fourvec Pmu_Hard_In = {0., 0., 0., 0.};
                fourvec Pmu_Hard_Out = {0., 0., 0., 0.};
                fourvec Pmu_Soft_Gain = {0., 0., 0., 0.};
                for (auto & p : plist) {
                    Pmu_Hard_In = Pmu_Hard_In + p.p;
                }
                while(med1.load_next()){
                    double current_hydro_clock = med1.get_tauL();
                    double hydro_dtau = med1.get_hydro_time_step();
                    //LOG_INFO << current_hydro_clock/fmc_to_GeV_m1 
                    //        << " [fm/c]\t" 
                    //       << " # of hard =" << plist.size();
                    // further divide hydro step into 10 transpor steps
                    int Ns = 1; 
                    double dtau = hydro_dtau/Ns, DeltaTau;
                    for (int i=0; i<Ns; ++i){
                        new_plist.clear();
                        for (auto & p : plist){     

                            if (p.Tf < 0.17) {
                                new_plist.push_back(p);
                                continue;       
                            }
                            if (p.x.tau() > current_hydro_clock+(i+1)*dtau){
                                // if the particle time is in the future 
                                // (not formed to the medium yet), put it back 
                                // in the list 
                                new_plist.push_back(p);
                                continue;
                            }
                            else{
                                DeltaTau = current_hydro_clock 
                                         + (i+1)*dtau - p.x.tau();
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
                             
                            Pmu_Soft_Gain = Pmu_Soft_Gain + ploss; 
                            current J; 
                            ploss = ploss.boost_to(0, 0, vzgrid);
                            J.p = ploss;
                            J.chetas = std::cosh(p.x.rap());
                            J.shetas = std::sinh(p.x.rap());
                            J.cs = std::sqrt(.3333);
                            clist.push_back(J);  
                        }
                        plist = new_plist;
                    }
                }
                for (auto & p : plist) {
                     if (std::abs(p.p.pseudorap())<3.)
                     Pmu_Hard_Out = Pmu_Hard_Out + p.p;
                }
                
                //LOG_INFO << "Hard initial " << Pmu_Hard_In;
                //LOG_INFO << "Hard final " << Pmu_Hard_Out;
                //LOG_INFO << "Soft deposite " << Pmu_Soft_Gain;
                std::stringstream fheader;
                fheader << args["output"].as<fs::path>().string() 
                        << processid;
                std::vector<double> Rs({.2,.4,.6,.8,1.});
                //clist.clear();
                FindJetTower(
                    plist, clist,
                    Rs, 10,
                    -3, 3,
                    fheader.str(), 
                    sigma_gen
                );
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



