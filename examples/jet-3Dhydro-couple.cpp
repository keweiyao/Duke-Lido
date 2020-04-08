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
        std::stringstream fheader;
        fheader << args["output"].as<fs::path>().string() 
                << processid;
        /// Initialize Lido in-medium transport
        initialize(table_mode,
            args["lido-setting"].as<fs::path>().string(),
            args["lido-table"].as<fs::path>().string()
        );
        /// Initialize a simple hadronizer
        JetDenseMediumHadronize Hadronizer;

        std::vector<double> pTbins({0,1,2,3,4,6,8,10,12,16,20,30,40,
               50,60,80,100,120,150,200,300,400,500,600,1000,2000});
        LeadingParton dNdpT(pTbins);
        // Scale to insert In medium transport
        double Q0 = args["Q0"].as<double>();
        std::vector<double> TriggerBin({
	5,10,15,20,25,30,40,50,60,80,100,
	110,120,130,140,150,160,170,180,200,
	220,240,260,280,300,
	350,400,450,500,550,600,650,700,
	750,800,900,1000,1200,1400,1600,
	1800,2000});
        for (int iBin = 0; iBin < TriggerBin.size()-1; iBin++){
            /// Initialize a pythia generator for each pT trigger bin
            PythiaGen pythiagen(
                args["pythia-setting"].as<fs::path>().string(),
                args["ic"].as<fs::path>().string(),
                TriggerBin[iBin],
                TriggerBin[iBin+1],
                args["eid"].as<int>(),
                Q0
            );
            
            LOG_INFO << pythiagen.sigma_gen() << " " 
                     << TriggerBin[iBin] << " "
                     << TriggerBin[iBin+1];
            for (int ie=0; ie<args["pythia-events"].as<int>(); ie++){
                std::vector<particle> plist, colorlist, hlist, thermal_list,
                                      new_plist, pOut_list;
                std::vector<current> clist;
                std::vector<HadronizeCurrent> slist;
                
                /// Initialzie a hydro reader
                Medium<3> med1(args["hydro"].as<fs::path>().string());
                // Initialize parton list from python
                pythiagen.Generate(plist);
                double sigma_gen = pythiagen.sigma_gen()
                                  / args["pythia-events"].as<int>();
                

                // freestream form t=0 to tau=tau0
                for (auto & p : plist){
                    p.tau_i=0.;
                    p.Tf = 0.161;
                    if (p.x.tau() < med1.get_tauH()){
                        p.freestream(
                            compute_realtime_to_propagate(
                                med1.get_tauH(), p.x, p.p
                            )
                        );
                    }
                }   

                // Energy-momentum checkbook
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
                            if (p.Tf < 0.16) {
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
                         
                            if (fs_size==-1){
                                // particle lost to the medium, but we
                                // track its color
                                colorlist.push_back(pOut_list[0]);
                            }
                            else {
                                for (auto & fp : pOut_list) {
                                    ploss = ploss - fp.p;
                                    new_plist.push_back(fp);
                                }
                            }               
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
                for (auto & p : colorlist) plist.push_back(p);
                Hadronizer.hadronize(plist, hlist, thermal_list, Q0, 1);
                for(auto & it : thermal_list){
                    //HadronizeCurrent J;
                    double vz = it.x.z()/it.x.t();
                    current J; 
                            J.p = it.p.boost_to(0, 0, vz)*(-1.);
                            J.chetas = std::cosh(it.x.rap());
                            J.shetas = std::sinh(it.x.rap());
                            J.cs = std::sqrt(.3333);
                            clist.push_back(J);  
                    /*double gamma = 1./std::sqrt(1.-vx*vx-vy*vy-vz*vz);
                    fourvec Umu{gamma, gamma*vx, gamma*vy, gamma*vz};
                    fourvec pmu{-it.p.t(), -it.p.x(), -it.p.y(), -it.p.z()};
                    J.UdotG = dot(Umu, pmu);
                    J.G = pmu;
                    J.gamma_perp = Umu.tau();
                    J.v_perp = std::sqrt(1.-1./(J.gamma_perp*J.gamma_perp));
                    J.etas = it.x.rap();
                    J.phiu = std::atan2(vy,vx);
                    slist.push_back(J);*/
                }

                std::vector<double> Rs({.2,.4,.6,.8, 1.});
                dNdpT.add_event(hlist, sigma_gen);
                /*FindJetTower(
                    plist, clist, slist,
                    Rs, 10,
                    -3, 3,
                    fheader.str(), 
                    sigma_gen
                );*/
            }
        }
        dNdpT.write(fheader.str());
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



