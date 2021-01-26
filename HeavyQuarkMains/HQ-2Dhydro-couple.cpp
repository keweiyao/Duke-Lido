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
#include "pythia_HQ_gen.h"
#include "Hadronize.h"

namespace po = boost::program_options;
namespace fs = boost::filesystem;


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
        /*if (!args.count("hydro")){
            throw po::required_option{"<hydro>"};
            return 1;
        }
        else{
            if (!fs::exists(args["hydro"].as<fs::path>())) {
                throw po::error{"<hydro> path does not exist"};
                return 1;
            }
        }*/

        std::vector<double> TriggerBin({
           2,3,4,6,8,10,12,14,16,18,20,22,24,26,28,30,
           35,40,45,50,55,60,65,70,75,80,85,90,100,
           110,120,140,160,180,200,
           240,280,320,360,400,500});

        /// Initialize Lido in-medium transport
	// Scale to insert In medium transport
        double Q0 = args["Q0"].as<double>();
        double muT = args["muT"].as<double>();
        double theta = args["theta"].as<double>();
	double cut = args["cut"].as<double>();
	double afix = args["afix"].as<double>();
        double Tf = args["Tf"].as<double>();
        std::vector<double> parameters{muT, afix, cut, theta};
        lido A(args["lido-setting"].as<fs::path>().string(), 
               args["lido-table"].as<fs::path>().string(), 
               parameters);
        A.set_frame(1); //Bjorken Frame
        HFHadronize VacShower;                
        /// Initialzie a hydro reader
        Medium<2> med1(args["hydro"].as<fs::path>().string());
        double mini_tau0 = med1.get_tauH();

        // Fill in all events
        LOG_INFO << "Events initialization, tau0 = " <<  mini_tau0;
        std::vector<particle> plist, dlist, flist;
        for (int iBin = 0; iBin < TriggerBin.size()-1; iBin++) {
            HQGenerator pythiagen(
                            args["pythia-setting"].as<fs::path>().string(),
                            args["ic"].as<fs::path>().string(),
                            TriggerBin[iBin],
                            TriggerBin[iBin+1],
                            args["eid"].as<int>(),
                            Q0
                            );
            dlist.clear();
            pythiagen.Generate(dlist,
                             args["pythia-events"].as<int>(),
                             4.0);
            plist.insert(plist.end(), dlist.begin(), dlist.end());
        }
        for (auto & p : plist){
            p.Tf = Tf+.001;
            if (p.x.x0() >= mini_tau0 ) continue;
            else {
                double dtau = std::max(mini_tau0, p.tau0);
                p.x.a[0] = dtau;
                p.x.a[1] += p.p.x()/p.p.t()*dtau;
                p.x.a[2] += p.p.y()/p.p.t()*dtau;
            }
        }

        while(med1.load_next()) {
            double current_hydro_clock = med1.get_tauL();
            double dtau = med1.get_hydro_time_step();
            LOG_INFO << "Hydro t = " 
                     << current_hydro_clock/5.076 << " fm/c";
            std::vector<particle> new_plist, pOut_list;
            for (auto & p : plist){     
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
		    p.p = put_on_shell(p.p, p.pid); 
                }
            }
        }
        // free some mem and boost to lab frame
	/*std::vector<particle> newlist, outlist;
	for (auto & p: plist){
     	    outlist.clear();
            VacShower.hadronize(p, outlist, Q0, p.Tf);
            for (auto ip:outlist) newlist.push_back(ip);
        }*/
 	for (auto & p : plist) {
            double tau = p.x.x0();
            double etas = p.x.x3();
            p.x.a[0] = tau*std::cosh(etas);
            p.x.a[3] = tau*std::sinh(etas);
            double vzcell = std::tanh(etas);
            double gamma = 1./std::sqrt(1.-vzcell*vzcell);
            double vz0cell = std::tanh(p.x0.x3());
            p.p = p.p.boost_back(0,0,vzcell);
            p.p0 = p.p0.boost_back(0,0,vz0cell);
            p.vcell[0] = p.vcell[0]/gamma/(1+vzcell*p.vcell[2]);
            p.vcell[1] = p.vcell[1]/gamma/(1+vzcell*p.vcell[2]);
            p.vcell[2] = (vzcell+p.vcell[2])/(1+vzcell*p.vcell[2]);

        }


        int processid = getpid();
        std::stringstream outputfilename1,outputfilename2;
        outputfilename1 << args["output"].as<fs::path>().string() 
                        << "c-quark-frzout.dat";
        outputfilename2 << args["output"].as<fs::path>().string()
                        << "b-quark-frzout.dat";
        output_oscar(plist, 4, outputfilename1.str());
        output_oscar(plist, 5, outputfilename2.str());
  
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



