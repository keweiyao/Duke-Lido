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

#include "Medium_Reader.h"
#include "simpleLogger.h"
#include "workflow.h"
#include "PGun.h"

namespace po = boost::program_options;
namespace fs = boost::filesystem;

void output_jet(std::string fname, std::vector<particle> plist){
    std::ofstream f(fname);
    for (auto & p : plist) f << p.pid << " " << p.p << std::endl;
    f.close();
}


int main(int argc, char* argv[]){
    using OptDesc = po::options_description;
    OptDesc options{};
    options.add_options()
          ("help", "show this help message and exit")
          ("nevents,n",
            po::value<int>()->value_name("INT")->default_value(1,"1"),
           "number of events")
          ("pid",
            po::value<int>()->value_name("INT")->default_value(21,"21"),
           "pid")
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
           ("muT",
           po::value<double>()->value_name("DOUBLE")->default_value(1.5,"1.5"),
           "mu_min/piT")
	   ("theta",
           po::value<double>()->value_name("DOUBLE")->default_value(4.,"4."),
           "Emin/T")
	   ("Tf",
           po::value<double>()->value_name("DOUBLE")->default_value(0.16,"0.16"),
           "Tf")
	   ("afix",
           po::value<double>()->value_name("DOUBLE")->default_value(-1.,"-1."),
           "fixed alpha_s, <0 for running alphas")
	   ("E0",
           po::value<double>()->value_name("DOUBLE")->default_value(100.,"100"),
           "E0")
	   ("cut",
           po::value<double>()->value_name("DOUBLE")->default_value(4.,"4."),
           "cut between diffusion and scattering, Qc^2 = cut*mD^2")
          ("ic,i",
            po::value<fs::path>()->value_name("PATH")->required(),
           "trento initial condition file")
          ("eid,j",
           po::value<int>()->value_name("INT")->default_value(0,"0"),
           "trento event id")
          ("hydro",
           po::value<fs::path>()->value_name("PATH")->required(),
           "hydro file")
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
            need_response_table = (fs::exists(args["response-table"].as<fs::path>())) ? false : true;
        }

        /// use process id to define filename
        int processid = getpid();
        
        /// Initialize Lido in-medium transport
	// Scale to insert In medium transport
        double muT = args["muT"].as<double>();
        double theta = args["theta"].as<double>();
	double cut = args["cut"].as<double>();
	double afix = args["afix"].as<double>();
        double Tf = args["Tf"].as<double>();

	initialize(table_mode,
            args["lido-setting"].as<fs::path>().string(),
            args["lido-table"].as<fs::path>().string(),
	    muT, afix, theta, cut
        );

        // Fill in all events
        std::vector<particle> plist;
        PGun pgungen(args["ic"].as<fs::path>().string());

        Medium<2> med1(args["hydro"].as<fs::path>().string());
        double mini_tau0 = med1.get_tauH();
        plist.clear();
        pgungen.Generate(args["pid"].as<int>(),
                            args["nevents"].as<int>(),
                            args["E0"].as<double>(),
                            plist);
        double E1 = 0.;
        for (auto & p : plist) E1 += p.p.t();
        E1 /= plist.size();    
        for (auto & p : plist){
                p.Tf = Tf+0.0001;
                p.origin = 0;
                if (p.x.tau() < mini_tau0 )
                    p.freestream(
			compute_realtime_to_propagate(
                                   mini_tau0 , p.x, p.p)
                        );
            }   
            std::vector<particle> pOut_list;
            
            while(med1.load_next()) {
              double current_hydro_clock = med1.get_tauL();
              for (int i=0;i<4;i++){
                double dtau = (i+1)*med1.get_hydro_time_step()/4.;
                for (auto & p : plist){     
                  if (p.Tf < Tf || std::abs(p.x.rap())>5. 			    || std::abs(p.p.rap()>5. )
		      ) {
                        continue;       
                  }
                  if (p.x.tau() > current_hydro_clock+dtau){
                    continue;
                  }
                  else{
                    double DeltaTau = current_hydro_clock + dtau - p.x.tau();
                    double T = 0.0, vx = 0.0, vy = 0.0, vz = 0.0;
                    med1.interpolate(p.x, T, vx, vy, vz);
                    int fs_size = update_particle_momentum_Lido(
                            DeltaTau, T, {vx, vy, vz}, 
                            p, pOut_list, Tf);          
                  }
                }
              }
            }
        double Efinal = 0., E2;
        for (auto & p : plist) {
             Efinal += p.p.t();
             E2 += p.p.t()*p.p.t();
        }
        Efinal /= plist.size();
        E2 /= plist.size();
        std::ofstream file("eloss.dat",std::ios_base::app);
        file << E1 << " " << E1-Efinal 
               << " " << std::sqrt((E2-Efinal*Efinal)/(plist.size()-1)) << std::endl;
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



