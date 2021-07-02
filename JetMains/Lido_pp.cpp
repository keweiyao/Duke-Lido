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

        std::vector<double> TriggerBin;
	if (args["jet"].as<bool>()){
           std::vector<double> a({10,20,30,40,50,60,70,80,90,100,110,
           120,140,160,180,200,
           220,240,260,280,300,400,500,600,800,1000,
           1500,2000,2500});


		TriggerBin = a;
	} else{
           std::vector<double> a({
           0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,
           35,40,45,50,55,60,70,80,90,100,110,120,130,
           150,170,190,210,230,250,300,350,400,450,500,
           600,700,800,1000,1500,2000,2500
                        });
	   TriggerBin = a;
	}

        std::vector<double> Rs({0.2, .4, 0.6, 0.8});

        std::vector<double> shaperbins({0., .05, .1, .15,  .2, .25, .3,
                         .35, .4, .45, .5,  .6, .7,  .8,
                          1., 1.5, 2.0, 2.5, 3.0});
	//std::vector<double> zbins({
	//		.00, .0125, .025, .05, .10, .15, .20, .25, .30, .35, .40, .45, 
	//		.50, .55, .60, .65, .70, .75, .80, .85, .90, .95, 1.0});
                std::vector<double> zbins({
                        0, 0.001, 0.002,
                         0.0035,.005,.0065,.0085, .01,
                        .012, .016, .02,.025, .03,.04,.05,
                        .07, .09, .120, .15, .20,
                        .25, .35, .45, .55, .65, .8,
                        1.
                        });

	std::vector<double> zpTbins({
         0.0, 0.15, 0.3, 0.5,   0.70626877,   0.99763116,   1.40919147,
         1.99053585,   2.81170663,   3.97164117,   5.61009227,
         7.92446596,  11.19360569,  15.8113883 ,  22.33417961,
         31.54786722,  44.56254691,  62.94627059,  88.9139705 ,
         125.59432158, 177.40669462, 250.59361681, 353.97289219,
         500.});

	LeadingParton HadronSample;
	JetStatistics JetSample(Rs, shaperbins, zbins, zpTbins);
        JetFinder jetfinder(300,300,3.);
        JetHFCorr jet_HF_corr(shaperbins);

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
                if (!pythiagen.Generate(e1.plist)) continue;
		e1.maxPT = pythiagen.maxPT();
                e1.sigma = pythiagen.sigma_gen()/args["pythia-events"].as<int>();
                e1.x0 = pythiagen.x0();		
                events.push_back(e1);
            }
        }

        // transform particles to lab frame
	for (auto & ie: events){
            for (auto & p : ie.plist) {
                p.p = p.p.boost_back(0,0,std::tanh(p.x.x3()));
            }
	    ie.hlist = ie.plist;
        }

        JetDenseMediumHadronize Hadronizer; 
	for (auto & ie : events){
            Hadronizer.hadronize(ie.plist, ie.hlist, ie.thermal_list,
                                 ie.Q0, 0.165, 1);
            if (args["jet"].as<bool>()){
                for(auto & it : ie.thermal_list){
                    double vz = std::tanh(it.x.x3());
                    current J;
                    J.p = it.p.boost_to(0, 0, vz)*(-1.);
                    J.etas = it.x.x3();
                    ie.clist.push_back(J);
                }
            }
        }
        LOG_INFO << "Jet finding in pp";
        /// use process id to define filename
        int processid = getpid();
        std::stringstream fheader;
        fheader << args["output"].as<fs::path>().string() 
                << "/lido";
        for (auto & ie : events){
            HadronSample.add_event(ie.hlist, ie.sigma);
            ie.clist.clear();
	    if (args["jet"].as<bool>()) {
                jetfinder.set_sigma(ie.sigma);
                jetfinder.MakeETower(
                     0.6, 0.165, args["pTtrack"].as<double>(),
                     ie.hlist, ie.clist, 10, true);
                jetfinder.FindJets(Rs, 10., -.9, .9, false);
                jetfinder.FindHF(ie.hlist);
                jetfinder.Frag(zbins, zpTbins);
		jetfinder.LabelFlavor();
                jetfinder.CalcJetshape(shaperbins);
	        JetSample.add_event(jetfinder.Jets, ie.sigma, ie.x0);
                jet_HF_corr.add_event(jetfinder.Jets, jetfinder.HFs, ie.sigma);
	    }
            ie.plist.clear();
        }
          LOG_INFO << "Jet finding done";    
        HadronSample.write(fheader.str());
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



