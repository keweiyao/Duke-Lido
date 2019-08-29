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

#include "simpleLogger.h"
#include "workflow.h"
#include "random.h"

namespace po = boost::program_options;
namespace fs = boost::filesystem;

int main(int argc, char* argv[]){
    using OptDesc = po::options_description;
    OptDesc options{};
    options.add_options()
          ("help", "show this help message and exit")
          ("nparticles",
            po::value<int>()->value_name("INT")->default_value(1000,"1000"),
           "number of hard partons")
          ("energy",
            po::value<double>()->value_name("FLOAT")->default_value(10.0,"10.0"),
           "fixed parton energy [GeV]")
          ("temp",
           po::value<double>()->value_name("FLOAT")->default_value(0.3,"0.3"),
           "medium temperature at t=t0 [GeV]")
          ("taui",
           po::value<double>()->value_name("FLOAT")->default_value(0.0,"0.0"),
           "medium starting time [fm/c]")
          ("tauf",
           po::value<double>()->value_name("FLOAT")->default_value(5.0,"5.0"),
           "medium stopping time [fm/c]")
          ("dtau",
           po::value<double>()->value_name("FLOAT")->default_value(0.1,"0.1"),
           "hydro time step dtau [fm/c]")
          ("nu",
           po::value<double>()->value_name("FLOAT")->default_value(0.5,"0.5"),
           "T = T0 * (tau0/tau) ^ (2-1/nu)")
          ("lido-setting,s", 
            po::value<fs::path>()->value_name("PATH")->required(),
           "Lido table setting file")
          ("lido-table,t", 
            po::value<fs::path>()->value_name("PATH")->required(),
           "Lido table path to file")
          ("mu,m", 
            po::value<double>()->value_name("FLOAT")->default_value(1.0,"1.0"),
            "medium scale paramtmer")
          ("afix,f", 
            po::value<double>()->value_name("FLOAT")->default_value(-1.0,"-1.0"),
            "fixed alphas value, -1 is running")
          ("k-factor,k", 
            po::value<double>()->value_name("FLOAT")->default_value(0.0,"0.0"),
            "K-factor for the delta-qhat")
          ("t-scale,a", 
            po::value<double>()->value_name("FLOAT")->default_value(1.0,"1.0"),
            "rescale the T-dependence")
          ("e-scale,b", 
            po::value<double>()->value_name("FLOAT")->default_value(1.0,"1.0"),
            "rescale the p-dependence")
          ("t-power,p", 
            po::value<double>()->value_name("FLOAT")->default_value(1.0,"1.0"),
            "T-dependence power")
          ("e-power,q", 
            po::value<double>()->value_name("FLOAT")->default_value(1.0,"1.0"),
            "p-dependence power")
          ("gamma,g", 
            po::value<double>()->value_name("FLOAT")->default_value(0.0,"0.0"),
            "kpara / kperp anisotropy parameter")
          ("qcut,c",
            po::value<double>()->value_name("FLOAT")->default_value(1.0,"1.0"), 
            "separation scale Q2 = qcut * mD2")  
          ("mass",
            po::value<double>()->value_name("FLOAT")->default_value(1.3,"1.3"), 
            "quark mass")  
          ("pid",
            po::value<int>()->value_name("INT")->default_value(4,"4"), 
            "particle pid, 4 or 21")  
          ("runid",
            po::value<int>()->value_name("INT")->default_value(0,"0"), 
            "running id")  
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
        // check pid
        if (!args.count("pid")){
            if ( (std::abs(args["pid"].as<int>()) != 4)
				&& (std::abs(args["pid"].as<int>()) != 21) ){
		        throw po::error{"pid only suppots 4, 21"};
		        return 1;
			}
        }

        // start
        double L = 5*5.026; // box size
        std::vector<particle> plist, new_plist, pOut_list;
		int pid = args["pid"].as<int>();
		int runid = args["runid"].as<int>();
		double E0 = args["energy"].as<double>();
		double M = args["mass"].as<double>();
        double pabs = std::sqrt(E0*E0-M*M);
        
		for (int i=0; i<args["nparticles"].as<int>(); i++){
		    double phi = Srandom::dist_phi(Srandom::gen);
		    double ctheta = Srandom::dist_costheta(Srandom::gen);
		    double stheta = std::sqrt(1.-ctheta*ctheta);
			fourvec p0{E0, pabs*stheta*std::cos(phi), pabs*stheta*std::sin(phi),pabs*ctheta };
            fourvec x0{0, 
                       L*Srandom::init_dis(Srandom::gen),
                       L*Srandom::init_dis(Srandom::gen),
                       L*Srandom::init_dis(Srandom::gen)
                      };
			particle entry;
			entry.p = p0; // 
			entry.p0 = p0; // 
			entry.mass = M; // mass
			entry.pid = Srandom::init_dis(Srandom::gen)>.5?pid:-pid; 
			entry.x0 = x0;
			entry.x = x0;
			entry.weight = 1.;
			entry.Tf = 0.0;
			entry.is_vac = false;
			entry.is_virtual = false;
			entry.vcell.resize(3);
			entry.vcell[0] = 0.; 
			entry.vcell[1] = 0.; 
			entry.vcell[2] = 0.; 
			plist.push_back(entry); 
		}

        /// Lido init
        initialize(table_mode,
                args["lido-setting"].as<fs::path>().string(),
                args["lido-table"].as<fs::path>().string(),
                args["mu"].as<double>(),
                args["afix"].as<double>(),
                args["k-factor"].as<double>(),
                args["t-scale"].as<double>(),
                args["e-scale"].as<double>(),
                args["t-power"].as<double>(),
                args["e-power"].as<double>(),
                args["gamma"].as<double>(),
                args["qcut"].as<double>(),
                0
                );

        /// Assign each quark a transverse position according to TRENTo Nbin output
        /// Freestream particle to the start of hydro
		double ti = args["taui"].as<double>() * 5.026; // convert to GeV^-1
		double tf = args["tauf"].as<double>() * 5.026; // convert to GeV^-1
		double dt = args["dtau"].as<double>() * 5.026; // convert to GeV^-1
		double T0 = args["temp"].as<double>();
		double nu = args["nu"].as<double>();
		int Nsteps = int((tf-ti)/dt)+1;
        for (auto & p : plist) p.freestream(ti);

        // initial E
        double Ei = mean_E(plist);

		std::ofstream f1("yield"+std::to_string(runid)+".dat");
		std::ofstream f2("tfinal_fdist"+std::to_string(runid)+".dat");
        for (int i=0; i<Nsteps; i++){ 
			double t = ti + dt*i;
            double T = T0 * std::pow(ti/t, 2./3.-1./nu/3.);
			if (i%100==0) LOG_INFO << t/5.026 << " [fm/c]\t" << T << " [GeV]" << " #=" << plist.size();

            //////// Free transport plus one (hard) particle evolution
            new_plist.clear();
            for (auto & p : plist){
                // for parton
                if (std::abs(p.pid) < 22)
                    OneBodyUpdate_Parton(dt, T, {0., 0., 0.}, p, pOut_list);  

                // for Hidden heavy flavor
                if (p.pid > 500)
                    int channel = OneBodyUpdate_Onium(dt, T, {0., 0., 0.}, p, pOut_list);

                for (auto & pp : pOut_list) new_plist.push_back(pp);
            }
            // update particle list, due to one (hard) particle evolution
            plist.clear();
            for (auto & p : new_plist) plist.push_back(p);

            // apply perodic boundary condition (due to free transport out of boundary)
            for (auto & p : plist) {       
                for (int k=1; k<4; k++){ 
                    if (p.x.a[k] > L) p.x.a[k] = p.x.a[k] - L;
                    if (p.x.a[k] < 0) p.x.a[k] = p.x.a[k] + L;
                }            
            }

            /////// Two (hard) particle evolution
            new_plist.clear();
            for (auto & p1 : plist){
                if (p1.pid == 5) {
                    auto xQ = p1.x;
                    for (auto & p2 : plist){
                        if (p2.pid == -5 && p2.is_virtual == false && p1.is_virtual == false) {
                            auto xQbar = p2.x;
                            double dist2 = 0.;
                            for (int k=1; k<4; k++){ 
                                double dx = std::abs(xQ.a[k]-xQbar.a[k]);
                                dist2 += std::pow(std::min(dx,L-dx), 2);
                            }
                            // for sufficiently close pair
                            if (dist2 < std::pow(1*5.026, 2)) {
                                for (int k=1; k<4; k++){
                                    if (std::abs(p1.x.a[k]-p2.x.a[k]) > L/2.){
                                        if (p1.x.a[k]>p2.x.a[k]) p2.x.a[k] += L;
                                        else p1.x.a[k] += L;
                                    }
                                }
                                int channel = TwoBodyUpdate_QQbar(dt, T, 
                                              {0., 0., 0.}, p1, p2, pOut_list);
                                if (channel != -1) {
                                    p1.is_virtual = true;
                                    p2.is_virtual = true;
                                    for (auto & pp : pOut_list) new_plist.push_back(pp);
                                    //LOG_INFO << "Recombined!";
                                }
                                else {
                                    p1.x = xQ;
                                    p2.x = xQbar;
                                }
                            }
                        }
                    }
                }
            }
            // update particle list, due to two (hard) particle evolution
            for(std::vector<particle>::iterator it=plist.begin();it!=plist.end();) {
                if (it->is_virtual) it = plist.erase(it);
                else it++;
            }
            for (auto & p : new_plist) plist.push_back(p);

            // apply perodic boundary condition (due to possible recombine)            
            for (auto & p : plist){
                for (int k=1; k<4; k++){ 
                    if (p.x.a[k] > L) p.x.a[k] = p.x.a[k] - L;
                    if (p.x.a[k] < 0) p.x.a[k] = p.x.a[k] + L;
                }
            }

            // output some statistics
            int N5 = 0,
                N553 = 0, 
                N100553 = 0, 
                N555 = 0;
                
            for (auto & p : plist) {
                if (p.pid == 555) N555++;
                if (p.pid == 553) N553++;
                if (p.pid == 100553) N100553++;
                if (std::abs(p.pid) == 5) N5 ++;
            }
            f1 << t/5.026 << "\t" << N5 << "\t"<< N553 
                           << "\t" << N100553 << "\t" << N555 << std::endl;
        }

        for (auto & p : plist)
           f2 << p.pid << " " << p.x << " " << p.p << std::endl;
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



