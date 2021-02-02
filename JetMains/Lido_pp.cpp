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
#include "lido.h"
#include "pythia_jet_gen.h"
#include "Hadronize.h"
#include "jet_finding.h"

namespace po = boost::program_options;
namespace fs = boost::filesystem;

void output_jet(std::string fname, std::vector<particle> plist){
    std::ofstream f(fname);
    f << "# pid, tau, x, y, etas, pT, phi, eta, M, Q0, color, acolor\n";
    for (auto & p : plist) f << p.pid << " " << p.x << " "
                             << p.p.xT() << " " << p.p.phi() << " "
                             << p.p.pseudorap() << " " << p.mass << " "
                             << p.Q0 << " "
                             << p.col << " " << p.acol
                             << std::endl;
    f.close();
}
struct event{
std::vector<particle> plist, thermal_list, hlist;
std::vector<current> clist;
double sigma, Q0, maxPT;
fourvec x0;
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
           ("output,o",
           po::value<fs::path>()->value_name("PATH")->default_value("./"),
           "output file prefix or folder")
	  ("jet", po::bool_switch(),
           "Turn on to do jet finding (takes time)")
	   ("pTtrack",
           po::value<double>()->value_name("DOUBLE")->default_value(.7,".7"),
           "minimum pT track in the jet shape reconstruction")
	   ("Q0,q",
           po::value<double>()->value_name("DOUBLE")->default_value(.5,".5"),
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

        std::vector<double> TriggerBin({
           0,2,4,6,8,10,12,14,16,18,20,24,28,32,
           36,40,45,50,60,70,80,90,100,110,
           120,140,160,180,200,250,300,350,400,500,
           600,700,800,900,1000,1200,1600,2000,2500});
	std::vector<double> Rs({.2,.3,.4, .6, .8,1.0});
        std::vector<double> ParticlepTbins({0,.25,.5,1.,1.5,2,3,4,5,6,7,8,10,
			12,14,16,18,20,22,24,28,32,36,40,45,50,
			55,60,65,70,80,90,100,
			110,120,140,160,180,200,220,240,260,300,
			350,400,500,600,800,1000});
        std::vector<double> jetpTbins({20,25,30,35,
			40,45,50,55,60,70,80,90,100,110,120,140,160,180,200,
			240,280,320,360,400,500,600,800,1000,
			1200,1400,1600,2000,2500});
        std::vector<double> HFpTbins({4,20,200});
        std::vector<double> HFETbins({2,6,10,20,40,100});
        std::vector<double> shapepTbins({60,80,100,120,2500});
        std::vector<double> shaperbins({0, .05, .1, .15,  .2, .25, .3,
                          .35, .4, .45, .5,  .6, .7,  .8,
                           1., 1.5, 2.0, 2.5, 3.0});
	std::vector<double> xJpTbins({100,126,158,178,200,224,251,282,316,398,562});
        std::vector<double> FragpTbins({40,60,70,100,126,158,200,251,316,398,600,800});

	std::vector<double> zbins({0, 0.001,0.002, 0.0035,.005,.0065,.0085,.011,.015,
                        .019,.025,.032,.042,.055,
                        .071, .092, .120,.157, .204, 
                        .266, .347, .452, .589, .767,
                        1.});
        std::vector<double> zpTbins({
         0.0, 0.15, 0.3, 0.5       ,   0.70626877,   0.99763116,   1.40919147,
         1.99053585,   2.81170663,   3.97164117,   5.61009227,
         7.92446596,  11.19360569,  15.8113883 ,  22.33417961,
         31.54786722,  44.56254691,  62.94627059,  88.9139705 ,
         125.59432158, 177.40669462, 250.59361681, 353.97289219,
         500.        	});

	LeadingParton dNdpT(ParticlepTbins);
	JetStatistics JetSample(jetpTbins, Rs, 
			 shapepTbins, shaperbins, 
			 FragpTbins, zbins, zpTbins,
			 xJpTbins);
        JetFinder jetfinder(300,300,3.);
        JetHFCorr jet_HF_corr(HFpTbins, shaperbins);

        std::vector<event> events;
        // Fill in all events
        double Q0 = args["Q0"].as<double>();
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
            for (int i=0; i<args["pythia-events"].as<int>(); i++){
                event e1;
		e1.Q0 = Q0;
                pythiagen.Generate(e1.plist);
		e1.maxPT = pythiagen.maxPT();
                e1.sigma = pythiagen.sigma_gen()/args["pythia-events"].as<int>();
	        LOG_INFO << pythiagen.sigma_gen() << " " <<e1.sigma;	
                e1.x0 = pythiagen.x0();		
                events.push_back(e1);
            }
        }

        // transform particles to lab frame
	for (auto & ie: events){
            for (auto & p : ie.plist) {
                p.p = p.p.boost_back(0,0,std::tanh(p.x.x3()));
            }
        }
      
        LOG_INFO << "Jet finding in pp";
        /// use process id to define filename
        int processid = getpid();
        std::stringstream fheader;
        fheader << args["output"].as<fs::path>().string() 
                << processid;
        for (auto & ie : events){
            dNdpT.add_event(ie.plist, ie.sigma, ie.maxPT);
            ie.clist.clear();
	    if (args["jet"].as<bool>()) {
                jetfinder.set_sigma(ie.sigma);
                jetfinder.MakeETower(
                     0.6, 0.15, args["pTtrack"].as<double>(),
                     ie.plist, ie.clist, 10, false);
                jetfinder.FindJets(Rs, 10., -3., 3., false);
                jetfinder.FindHF(ie.plist);
                jetfinder.Frag(zbins, zpTbins);
		jetfinder.LabelFlavor();
                jetfinder.CalcJetshape(shaperbins);
	        JetSample.add_event(jetfinder.Jets, ie.sigma, ie.x0);
                jet_HF_corr.add_event(jetfinder.Jets, jetfinder.HFs, ie.sigma);
	    }
            ie.plist.clear();
        }
          LOG_INFO << "Jet finding idone";    
        dNdpT.write(fheader.str());
	if (args["jet"].as<bool>()){
           JetSample.write(fheader.str());
           jet_HF_corr.write(fheader.str());
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



