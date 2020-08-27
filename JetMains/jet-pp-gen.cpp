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
           ("output,o",
           po::value<fs::path>()->value_name("PATH")->default_value("./"),
           "output file prefix or folder")
          ("Q0,q",
           po::value<double>()->value_name("DOUBLE")->default_value(.4,".4"),
           "Scale [GeV] to insert in-medium transport")
	  ("jet", po::bool_switch(),
           "Turn on to do jet finding (takes time)")
	  ("pTtrack",
           po::value<double>()->value_name("DOUBLE")->default_value(.7,".7"),
           "minimum pT track in the jet shape reconstruction")

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
/*
        /// all kinds of bins and cuts
        // For RHIC 200 GeV
        std::vector<double> TriggerBin({
         2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,18,20,22,24,28,32,36,40,45,50,55,60,70,80,90,100});
        std::vector<double> Rs({.2,.3,.4});
        std::vector<double> ParticlepTbins({0,1,2,3,4,5,6,8,10,12,14,16,18,20,22,24,26,30,40,60,100});
        std::vector<double> jetpTbins({3,4,5,6,7,10,12,14,16,18,20,24,28,32,36,40,50,60,100});
        std::vector<double> HFpTbins({2,6,10,20,40,100});
        std::vector<double> HFETbins({2,6,10,20,40,100});
        std::vector<double> shapepTbins({20,30,40,60,80,120,2000});
        std::vector<double> shaperbins({0, .05, .1, .15,  .2, .25, .3,
                          .35, .4, .45, .5,  .6, .7,  .8,
                           1., 1.5, 2.0, 2.5, 3.0});
        std::vector<double> xJpTbins({8,12,16,20,30,40,60});
        std::vector<double> FragpTbins({10,20,30,40});
        std::vector<double> zbins({.005,.0065,.0085,.011,.015,
                        .019,.025,.032,.042,.055,
                        .071, .092, .120,.157, .204,
                        .266, .347, .452, .589, .767,
                        1.});

*/	
	
	// For 5.02 TeV
        std::vector<double> TriggerBin({
         2,4,6,8,10,12,14,16,20,
         24,28,32,36,40,50,60,70,80,90,100,
         110,120,130,140,150,160,180,200,240,280,320,360,400,500,
         600,700,800,1000,1200,1500,2000,2500});
        std::vector<double> Rs({0.2,0.3,.4,.6,.8,1.0});
        std::vector<double> ParticlepTbins({0,1,2,3,4,5,6,8,10,12,14,16,20,
                        24,28,32,40,50,60,70,80,90,100,110,
                        120,130,140,150,160,180,200,250,300,350,400,600,800,1000});
        std::vector<double> jetpTbins({4,6,8,10,12,15,20,25,30,
                        40,50,60,70,80,100,120,140,160,180,200,
                        240,280,320,360,400,500,600,800,1000,
                        1200,1400,1600,2000,2500});
        std::vector<double> HFpTbins({4,20,200});
        std::vector<double> HFETbins({2,6,10,20,40,100});
        std::vector<double> shapepTbins({60,80,100,120,2000});
        std::vector<double> shaperbins({0, .05, .1, .15,  .2, .25, .3,
                          .35, .4, .45, .5,  .6, .7,  .8,
                           1., 1.5, 2.0, 2.5, 3.0});
        std::vector<double> xJpTbins({100,126,158,178,200,224,251,282,316,398,562});
        std::vector<double> FragpTbins({100,126,158,200,251,316,398,600,800});
        std::vector<double> zbins({.005,.0065,.0085,.011,.015,
                        .019,.025,.032,.042,.055,
                        .071, .092, .120,.157, .204, 
                        .266, .347, .452, .589, .767,
                        1.});

        LeadingParton dNdpT(ParticlepTbins);
        JetStatistics JetSample(jetpTbins, Rs, 
			shapepTbins, shaperbins, 
			FragpTbins, zbins,
			xJpTbins);
        JetHFCorr jet_HF_corr(HFpTbins, shaperbins);
        HFETCorr  HF_ET_corr(HFETbins, shaperbins);
        /// Initialize jet finder with medium response
        JetFinder jetfinder(150,150,3.);


        // Scale to insert In medium transport
        double Q0 = args["Q0"].as<double>();
        /// use process id to define filename
        int processid = getpid();
        std::stringstream fheader;
        fheader << args["output"].as<fs::path>().string() 
                << processid;
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
            //LOG_INFO << pythiagen.sigma_gen() << " "
            //         << TriggerBin[iBin] << " "
            //         << TriggerBin[iBin];
            for (int ie=0; ie<args["pythia-events"].as<int>(); ie++){
                std::vector<particle> plist;
                std::vector<current> clist;
                std::vector<HadronizeCurrent> slist;
                // Initialize parton list from python
                pythiagen.Generate(plist);
                double sigma_gen = pythiagen.sigma_gen()
                                  / args["pythia-events"].as<int>();
                
                dNdpT.add_event(plist, sigma_gen, pythiagen.maxPT());
	        if (args["jet"].as<bool>()) {
                    jetfinder.set_sigma(sigma_gen);
                    jetfinder.MakeETower(
                         0.6, 0.165, args["pTtrack"].as<double>(),
                         plist, clist, slist, 10);
                    jetfinder.FindJets(Rs, 10., -3., 3.);
                    jetfinder.FindHF(plist);
                    jetfinder.CorrHFET(shaperbins);
                    jetfinder.LabelFlavor();
                    jetfinder.CalcJetshape(shaperbins);
		    jetfinder.Frag(zbins);
	            JetSample.add_event(jetfinder.Jets, sigma_gen, pythiagen.x0());
                    jet_HF_corr.add_event(jetfinder.Jets, jetfinder.HFs,
                                          sigma_gen);
                    HF_ET_corr.add_event(jetfinder.HFaxis, sigma_gen);
	        }
            }
        }
        dNdpT.write(fheader.str());
	if (args["jet"].as<bool>()){
	    JetSample.write(fheader.str());
	    jet_HF_corr.write(fheader.str());
	    HF_ET_corr.write(fheader.str());
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



