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
          ("pthat-low,l",
            po::value<int>()->value_name("INT")->default_value(100,"100"),
           "pThatMin")
          ("pthat-high,h",
            po::value<int>()->value_name("INT")->default_value(120,"120"),
           "pThatMax")
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
        std::ofstream fb("balance.txt");
        /// Initialize pythia generator
        PythiaGen pythiagen(
                args["pythia-setting"].as<fs::path>().string(),
                args["ic"].as<fs::path>().string(),
                args["pthat-low"].as<int>(),
                args["pthat-high"].as<int>(),
                args["eid"].as<int>()
        );
        /// Initialize Lido in-medium transport
        LOG_INFO<<"init lido";
        initialize(table_mode,
                args["lido-setting"].as<fs::path>().string(),
                args["lido-table"].as<fs::path>().string()
                );
        /// Initialize a simple hadronizer
        JetDenseMediumHadronize Hadronizer;
        
        std::ostringstream outname1, outname2;
        outname1 << "I-" << args["pthat-low"].as<int>()
                 << "-" << args["pthat-high"].as<int>();
        outname2 << "F-" << args["pthat-low"].as<int>()
                 << "-" << args["pthat-high"].as<int>();
        for (int ie=0; ie<args["pythia-events"].as<int>(); ie++){
            std::vector<particle> plist, hlist, thermal_list,
                                  new_plist, pOut_list;
            std::vector<current> clist;
            // Initialize parton list from python
            pythiagen.Generate(plist, args["heavy"].as<int>());
            double sigma_gen = plist[0].weight;
          /* for (int ir=0; ir<5; ir++){
                std::ostringstream ss;
                ss << outname1.str() << "-Raa-R" << (ir+1) << ".dat";
                FindJet(plist, clist,
                 0.2*(ir+1), 60, -2.8, 2.8,
                 ss.str(), sigma_gen);
            }continue;  

            std::ostringstream sss;
                sss << outname1.str() << "-R0d4" << ".dat";
                JetShape(plist, clist,
                 0.4, 80, -2, 2,
                 sss.str(), sigma_gen);
            continue;*/
            /// Initialzie a hydro reader
            Medium<2> med1(args["hydro"].as<fs::path>().string());

            // freestream form t=0 to tau=tau0
            for (auto & p : plist){
                p.Tf = 0.155;
                p.origin=-1;
                if (p.x.tau() < med1.get_tauH()){
                    p.freestream(compute_realtime_to_propagate(
                                 med1.get_tauH(), p.x, p.p));
                }else{
                    LOG_INFO << "need to wait";
                }
            }   

            // Energy-momentum checkbook
            fourvec Pmu_Hard_In = {0., 0., 0., 0.};
            fourvec Pmu_Hard_Out = {0., 0., 0., 0.};
            fourvec Pmu_Soft_Gain = {0., 0., 0., 0.};
            std::ofstream fsoft("soft.txt");
            for (auto & p : plist) {
                if (std::abs(p.p.pseudorap())<2)
                    Pmu_Hard_In = Pmu_Hard_In + p.p;
            }
            while(med1.load_next()){
                double current_hydro_clock = med1.get_tauL();
                double hydro_dtau = med1.get_hydro_time_step();
                LOG_INFO << current_hydro_clock/fmc_to_GeV_m1 
                         << " [fm/c]\t" 
                         << " # of hard =" << plist.size();
                // further divide hydro step into 10 transpor steps
                int Ns = 2; 
                double dtau = hydro_dtau/Ns, DeltaTau;
                for (int i=0; i<Ns; ++i){
                    new_plist.clear();
                    for (auto & p : plist){     
                        if (p.Tf < 0.15) {
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
                        double vzgrid = p.x.z()/p.x.t();
                        med1.interpolate(p.x, T, vx, vy, vz);

                        fourvec ploss = p.p;
                        int fs_size = update_particle_momentum_Lido(
                                  DeltaTau, T, {vx, vy, vz}, p, pOut_list);
                        
                        for (auto & fp : pOut_list) {
                            ploss = ploss - fp.p;
                            new_plist.push_back(fp);
                        }
                       if (std::abs(ploss.pseudorap())<2)            
                           Pmu_Soft_Gain = Pmu_Soft_Gain + ploss;
                   
                    current J; 
                    vx=0; vy=0; vz=vzgrid;
                    ploss = ploss.boost_to(vx, vy, vz); 
                    J.p = ploss;
                    J.x = p.x;
                    J.v[0] = vx; J.v[1] = vy; J.v[2] = vz;
                    clist.push_back(J);           
                   
                    }
                    plist = new_plist;
                }
            }
            for (auto & p : plist) {
                if (std::abs(p.p.pseudorap())<2)
                Pmu_Hard_Out = Pmu_Hard_Out + p.p;
            }
            Hadronizer.hadronize(plist, hlist, thermal_list);
          for(auto & it : thermal_list){
                current J;
                fourvec pmu{-it.p.t(), -it.p.x(), -it.p.y(), -it.p.z()};
                J.p = pmu.boost_to(it.vcell[0], it.vcell[1], it.vcell[2]);
                J.x = it.x;
                J.v[0] = it.vcell[0];
                J.v[1] = it.vcell[1];
                J.v[2] = it.vcell[2];
                clist.push_back(J);
            }
              for (int ir=0; ir<5; ir++){
                std::ostringstream ss;
                ss << outname2.str() << "-Raa-R" << (ir+1) << ".dat";
                FindJet(plist, clist,
                 0.2*(ir+1), 60, -2.8, 2.8,
                 ss.str(), sigma_gen);
            } /*
             std::ostringstream ss;
                ss << outname2.str() << "-R0d4" << ".dat";
                JetShape(hlist, clist,
                 0.4, 80, -2, 2,
                 ss.str(), sigma_gen);
            */
            //
            LOG_INFO << "Hard initial " << Pmu_Hard_In;
            LOG_INFO << "Hard final " << Pmu_Hard_Out;
            LOG_INFO << "Soft deposite " << Pmu_Soft_Gain;

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



