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

struct event{
std::vector<particle> plist, colorlist, thermal_list, hlist;
std::vector<current> clist;
std::vector<HadronizeCurrent> slist;
double sigma;
};

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

        /// all kinds of bins and cuts
	std::vector<double> TriggerBin({
         1,3,5,10,15,20,30,40,60,
         80,100,120,
         140,160,180,200,240,280,
         320,360,400,450,500,550,
         600,700,800,900,
         1000,1200,1400,1800,2500});

        std::vector<double> Rs({.2,.4,.6,.8, 1.});

        std::vector<double> ParticlepTbins({0,1,2,3,4,6,8,10,12,16,20,30,40,
               50,60,80,100,120,150,200,300,400,500,600,1000,2000});
        std::vector<double> jetpTbins({10,15,20,30,40,
               50,60,80,120,160,200,300,400,500,600,800,1000,1400,1800});
	std::vector<double> shapepTbins({20,30,40,60,80,120,2000});
        std::vector<double> shaperbins({0, .05, .1, .15,  .2, .25, .3,
                          .35, .4, .45, .5,  .6, .7,  .8,
                           1., 1.5, 2.0, 2.5, 3.0});

	LeadingParton dNdpT(ParticlepTbins);
	JetStatistics JetSample(jetpTbins, Rs, shapepTbins, shaperbins);
        // Scale to insert In medium transport
        double Q0 = args["Q0"].as<double>();
                
        /// Initialzie a hydro reader
        Medium<2> med1(args["hydro"].as<fs::path>().string());
        std::vector<event> events;
        // Fill in all events
        LOG_INFO << "Events initialization";
        for (int iBin = 0; iBin < TriggerBin.size()-1; iBin++) {    
            /// Initialize a pythia generator for each pT trigger bin
            PythiaGen pythiagen(
                args["pythia-setting"].as<fs::path>().string(),
                args["ic"].as<fs::path>().string(),
                TriggerBin[iBin],
                TriggerBin[iBin+1],
                args["eid"].as<int>(),
                Q0
            );
            LOG_INFO << " Generating " 
                     << TriggerBin[iBin] << " < PT_hat <"
                     << TriggerBin[iBin+1] << " GeV";
            for (int i=0; i<args["pythia-events"].as<int>(); i++){
                event e1;
                pythiagen.Generate(e1.plist);
                e1.sigma = pythiagen.sigma_gen()
                           /args["pythia-events"].as<int>();            
                // freestream form t=0 to tau=tau0
                for (auto & p : e1.plist){
                    p.Tf = 0.161;
                    if (p.x.tau() < med1.get_tauH())
                        p.freestream(compute_realtime_to_propagate(
                                   med1.get_tauH(), p.x, p.p)
                        );
                }   
                events.push_back(e1);
            }
        }
        LOG_INFO << "Start evolution of " << events.size() << " hard events";
        while(med1.load_next()) {
            double current_hydro_clock = med1.get_tauL();
            double dtau = med1.get_hydro_time_step();
            LOG_INFO << "Hydro t = " << current_hydro_clock/5.076 << " fm/c";
            for (auto & ie : events){
                std::vector<particle> new_plist, pOut_list;
                for (auto & p : ie.plist){     
                    if (p.Tf < 0.15 || std::abs(p.p.rap())>5.) {
                        new_plist.push_back(p);
                        continue;       
                    }
               
                    if (p.x.tau() > current_hydro_clock+dtau){
                        // if the particle time is in the future 
                        // (not formed to the medium yet), put it back 
                        // in the list 
                        new_plist.push_back(p);
                        continue;
                    }
                    else{
                        double DeltaTau = 
                               current_hydro_clock + dtau - p.x.tau();
                        // get hydro information
                        double T = 0.0, vx = 0.0, vy = 0.0, vz = 0.0;
                        med1.interpolate(p.x, T, vx, vy, vz);
                        fourvec ploss = p.p;
                        int fs_size = update_particle_momentum_Lido(
                                DeltaTau, T, {vx, vy, vz}, p, pOut_list);     
                        if (fs_size==-1){
                            // particle lost to the medium, but we
                            // track its color
                            //ie.colorlist.push_back(pOut_list[0]);
			    new_plist.push_back(pOut_list[0]);
			    current J;
                            ploss = ploss - pOut_list[0].p;
			    ploss = ploss.boost_to(0, 0, p.x.z()/p.x.t());
			    J.p = ploss;
                            J.chetas = std::cosh(p.x.rap());
                            J.shetas = std::sinh(p.x.rap());
                            J.cs = std::sqrt(.3333);
                            ie.clist.push_back(J);
                        }
                        else {
                            for (auto & fp : pOut_list) {
                                ploss = ploss - fp.p;
                                new_plist.push_back(fp);
                            }                           
                            current J; 
                            ploss = ploss.boost_to(0, 0, p.x.z()/p.x.t());
                            J.p = ploss;
                            J.chetas = std::cosh(p.x.rap());
                            J.shetas = std::sinh(p.x.rap());
                            J.cs = std::sqrt(.3333);
                            ie.clist.push_back(J);  
                        }
                    }
                }
                ie.plist = new_plist;
            }
        }    

        // Hadronization
        // put back lost particles with their color
        LOG_INFO << "Hadronization";
        for (auto & ie : events){
            for (auto & p : ie.colorlist) ie.plist.push_back(p);
            Hadronizer.hadronize(ie.plist, ie.hlist, ie.thermal_list, Q0, 1);
            for(auto & it : ie.thermal_list){
                double vz = it.x.z()/it.x.t();
                current J; 
                J.p = it.p.boost_to(0, 0, vz)*(-1.);
                J.chetas = std::cosh(it.x.rap());
                J.shetas = std::sinh(it.x.rap());
                J.cs = std::sqrt(.3333);
                ie.clist.push_back(J);  
            }
         }
        LOG_INFO << "Jet finding, w/ medium excitation";
        for (auto & ie : events){
            dNdpT.add_event(ie.hlist, ie.sigma);
            auto jets = FindJetTower(
                 ie.hlist, ie.clist, ie.slist, 
	         Rs, shaperbins, 10, -3, 3, ie.sigma);
	    JetSample.add_event(jets, ie.sigma);
        }
        dNdpT.write(fheader.str());
	JetSample.write(fheader.str());
	
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



