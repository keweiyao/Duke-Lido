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
          ("hydro",
           po::value<fs::path>()->value_name("PATH")->required(),
           "hydro file")
          ("lido-setting,s", 
            po::value<fs::path>()->value_name("PATH")->required(),
           "Lido table setting file")
          ("lido-table,t", 
            po::value<fs::path>()->value_name("PATH")->required(),
           "Lido table path to file")  
          ("response-table,r", 
            po::value<fs::path>()->value_name("PATH")->required(),
           "response table path to file")  
           ("output,o",
           po::value<fs::path>()->value_name("PATH")->default_value("./"),
           "output file prefix or folder")
	  ("jet", po::bool_switch(),
           "Turn on to do jet finding (takes time)")
	   ("pTtrack",
           po::value<double>()->value_name("DOUBLE")->default_value(.7,".7"),
           "minimum pT track in the jet shape reconstruction")
           ("muT",
           po::value<double>()->value_name("DOUBLE")->default_value(1.5,"1.5"),
           "mu_min/piT")
	   ("Q0,q",
           po::value<double>()->value_name("DOUBLE")->default_value(.4,".4"),
           "Scale [GeV] to insert in-medium transport")
	   ("theta",
           po::value<double>()->value_name("DOUBLE")->default_value(4.,"4."),
           "Emin/T")
	   ("afix",
           po::value<double>()->value_name("DOUBLE")->default_value(-1.,"-1."),
           "fixed alpha_s, <0 for running alphas")
	   ("cut",
           po::value<double>()->value_name("DOUBLE")->default_value(4.,"4."),
           "cut between diffusion and scattering, Qc^2 = cut*mD^2")
           ("Tf",
           po::value<double>()->value_name("DOUBLE")->default_value(0.17,"0.17"),
           "Transport stopping temperature, Tf")
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
        // check lido-response table
        bool need_response_table;
        if (!args.count("response-table")){
            throw po::required_option{"<response-table>"};
            return 1;
        }
        else{
            need_response_table = (fs::exists(args["response-table"].as<fs::path>())) ? 
                         false : true;
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

        std::vector<double> TriggerBin({
           2,4,6,8,10,12,16,20,
           30,40,50,60,70,80,100,
           120,140,160,180,200,250,300,350,400,500,
           600,700,800,900,1000,1200,1600,2000,2500});
	std::vector<double> Rs({.2,.3,.4});
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
        JetFinder jetfinder(300,300,3., need_response_table, args["response-table"].as<fs::path>().string());
        JetHFCorr jet_HF_corr(HFpTbins, shaperbins);

        /// Initialize Lido in-medium transport
	// Scale to insert In medium transport
        double Q0 = args["Q0"].as<double>();
        double muT = args["muT"].as<double>();
        double theta = args["theta"].as<double>();
	double cut = args["cut"].as<double>();
	double afix = args["afix"].as<double>();
        double Tf = args["Tf"].as<double>();
        std::vector<double> parameters{muT, afix, cut, theta};
        JetDenseMediumHadronize Hadronizer;
        lido A(args["lido-setting"].as<fs::path>().string(), 
               args["lido-table"].as<fs::path>().string(), 
               parameters);
        A.set_frame(1); //Bjorken Frame
                
        /// Initialzie a hydro reader
        Medium<2> med1(args["hydro"].as<fs::path>().string());
        double mini_tau0 = med1.get_tauH();
        std::vector<event> events;
        // Fill in all events
        LOG_INFO << "Events initialization, tau0 = " <<  mini_tau0;
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
                e1.sigma = pythiagen.sigma_gen()
                           /args["pythia-events"].as<int>();  
                e1.x0 = pythiagen.x0();		
                // freestream form t=0 to tau=tau0
                // move partciles below tau0 to tau0
                for (auto & p : e1.plist){
                    p.Tf = Tf+.001;
                    if (p.x.x0() >= mini_tau0 ) continue;
                    else {
                        double dtau = std::max(mini_tau0, p.tau0);
                        p.x.a[0] = dtau;
                        p.x.a[1] += p.p.x()/p.p.t()*dtau;
                        p.x.a[2] += p.p.y()/p.p.t()*dtau;
                    }
                }   
                events.push_back(e1);
            }
        }

        LOG_INFO << "Start evolution of " << events.size() << " hard events";
        while(med1.load_next()) {
            double current_hydro_clock = med1.get_tauL();
            double dtau = med1.get_hydro_time_step();
            LOG_INFO << "Hydro t = " 
                     << current_hydro_clock/5.076 << " fm/c";
            for (auto & ie : events){
                std::vector<particle> new_plist, pOut_list;
                for (auto & p : ie.plist){     
                    if ( std::abs(p.x.x3())>5. || p.Tf<Tf ){
                        // skip particles at large space-time rapidity
                        new_plist.push_back(p);
                        continue;       
                    }
                    if (p.x.x0() > current_hydro_clock+dtau){
                        // skip particles in the future
                        new_plist.push_back(p);
                        continue;
                    }
                    else{
                        double DeltaTau = current_hydro_clock 
                                        + dtau - p.x.x0();
                        fourvec ploss = p.p;
                        double T = 0.0, vx = 0.0, vy = 0.0, vz = 0.0;
                        med1.interpolate(p.x, T, vx, vy, vz);
                        pOut_list.clear();
                        A.update_single_particle(DeltaTau, 
                                                 T, {vx, vy, vz}, 
                                                 p, pOut_list
                                                 );  
                        for (auto & fp : pOut_list) {
                            // compute energy momentum loss of hard partons
                            ploss = ploss - fp.p;
                            new_plist.push_back(fp);
                        }
                        if (args["jet"].as<bool>()){
                            current J; 
                            J.p = ploss;
                            J.etas = p.x.x3();
                            ie.clist.push_back(J);  
                        }
                    }
                }
                ie.plist = new_plist;
	    }
        }
        // free some mem
	for (auto & ie: events)
            for (auto & p : ie.plist) 
                p.radlist.clear();

        for (auto & ie : events){ 
            Hadronizer.hadronize(ie.plist, ie.hlist, ie.thermal_list,
                                 ie.Q0, Tf, 1);
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
        LOG_INFO << "Jet finding, w/ medium excitation";
        /// use process id to define filename
        int processid = getpid();
        std::stringstream fheader;
        fheader << args["output"].as<fs::path>().string() 
                << processid;
        for (auto & ie : events){
            dNdpT.add_event(ie.hlist, ie.sigma, ie.maxPT);
	    if (args["jet"].as<bool>()) {
                jetfinder.set_sigma(ie.sigma);
                jetfinder.MakeETower(
                     0.6, Tf, args["pTtrack"].as<double>(),
                     ie.hlist, ie.clist, 10, false);
                jetfinder.FindJets(Rs, 5., -3., 3., false);
                jetfinder.FindHF(ie.hlist);
                jetfinder.Frag(zbins, zpTbins);
		jetfinder.LabelFlavor();
                jetfinder.CalcJetshape(shaperbins);
	        JetSample.add_event(jetfinder.Jets, ie.sigma, ie.x0);
                jet_HF_corr.add_event(jetfinder.Jets, jetfinder.HFs, ie.sigma);
	    }
        }
        dNdpT.write(fheader.str());
	if (args["jet"].as<bool>()){
            LOG_INFO << fheader.str();
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



