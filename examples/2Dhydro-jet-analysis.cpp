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
//#include "jet_finding.h"

#include "JetCrossSectionAnalyzer.h"

namespace po = boost::program_options;
namespace fs = boost::filesystem;

void output_jet(std::string fname, std::vector<particle> plist){
    std::ofstream f(fname);
    for (auto & p : plist) f << p.pid << " " << p.p << " " << p.x0 << " " << p.x << std::endl;
    f.close();
}

int main(int argc, char* argv[]){

    std::vector<double> pTBin({39, 50, 63, 79, 100, 125, 158, 199, 251, 316});
	JetCrossSectionAnalyzer analyzer=JetCrossSectionAnalyzer(pTBin, -1, 0.4, 2, 0, 2.1);

    std::vector<double> pTHatBin({10,20,50,70,100,150,200,300,400,500,700,1000});
	size_t nBin = pTHatBin.size() - 1;
	
    for (size_t iBin = 0; iBin < nBin; iBin++)
	{

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
        
        /// Initialize Lido in-medium transport
        LOG_INFO<<"init lido";
        initialize(table_mode,
                args["lido-setting"].as<fs::path>().string(),
                args["lido-table"].as<fs::path>().string()
                );
        /// Initialize a simple hadronizer
        JetDenseMediumHadronize Hadronizer;


        std::ofstream fb("balance.txt");
        /// Initialize pythia generator
        PythiaGen pythiagen(
                args["pythia-setting"].as<fs::path>().string(),
                args["ic"].as<fs::path>().string(),
                pTHatBin[iBin],
                pTHatBin[iBin+1],
                args["eid"].as<int>()
        );

        for (int ie=0; ie<args["pythia-events"].as<int>(); ie++){
            std::vector<particle> plist, hlist, thermal_list,
                                  new_plist, pOut_list;
            std::vector<fourvec> clist;
            // Initialize parton list from python
            pythiagen.Generate(plist, args["heavy"].as<int>());
            double sigma_gen = plist[0].weight;

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
                if (std::abs(p.p.rap())<2)
                    Pmu_Hard_In = Pmu_Hard_In + p.p;
            }
            while(med1.load_next()){
                double current_hydro_clock = med1.get_tauL();
                double hydro_dtau = med1.get_hydro_time_step();
                //LOG_INFO << current_hydro_clock/fmc_to_GeV_m1 
                //         << " [fm/c]\t" 
                //         << " # of hard =" << plist.size();
                // further divide hydro step into 10 transpor steps
                int Ns = 10; 
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
                       if (std::abs(ploss.rap())<2)            
                       Pmu_Soft_Gain = Pmu_Soft_Gain + ploss;
                       //if (ploss.pabs() > .3)
                       //clist.push_back(ploss);
                    }
                    plist = new_plist;
                }
            }

            std::vector<fastjet::PseudoJet> fjInputs;
            for (auto & p : plist) {
                if (std::abs(p.p.rap())<2)
                Pmu_Hard_Out = Pmu_Hard_Out + p.p;
                fjInputs.push_back(fastjet::PseudoJet(p.p.x(), p.p.y(), p.p.z(), p.p.t()));
            }

            double sigmapb_weight = sigma_gen * 1.0e9/args["pythia-events"].as<int>();

			analyzer.doJetFinding(fjInputs, sigmapb_weight);
            
            //Hadronizer.hadronize(plist, hlist, thermal_list);
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
    }

    analyzer.outputResults();

    return 0;
}



