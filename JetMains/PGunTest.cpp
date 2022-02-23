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
#include "PGunWithShower.h"
#include "Hadronize.h"
#include "jet_finding.h"

namespace po = boost::program_options;
namespace fs = boost::filesystem;

void output_jet(std::string fname, std::vector<particle> plist){
    std::ofstream f(fname);
    f << "# pid, tau, x, y, etas, pT, phi, eta, M\n";
    for (auto & p : plist) f << p.pid << " "
                             << p.p.xT() << " " << p.p.phi() << " "
                             << p.p.pseudorap()  << std::endl;
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
       ("Q0,q",
           po::value<double>()->value_name("DOUBLE")->default_value(.5,".5"),
           "Scale [GeV] to insert in-medium transport")
       ("theta",
           po::value<double>()->value_name("DOUBLE")->default_value(4.,"4."),
           "Emin/T")
       ("Tf",
           po::value<double>()->value_name("DOUBLE")->default_value(0.16,"0.16"),
           "Tf")
       ("pTmin",
           po::value<double>()->value_name("DOUBLE")->default_value(0.7,"0.7"),
           "pTmin")
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
          ("pid",
           po::value<int>()->value_name("INT")->default_value(1,"1"),
           "hard id")
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
        double Q0 = args["Q0"].as<double>();
        double muT = args["muT"].as<double>();
        double theta = args["theta"].as<double>();
        double cut = args["cut"].as<double>();
        double afix = args["afix"].as<double>();
        double Tf = args["Tf"].as<double>();
        int PID = args["pid"].as<int>();
        std::vector<double> parameters{muT, afix, cut, theta};
        lido A(args["lido-setting"].as<fs::path>().string(), 
               args["lido-table"].as<fs::path>().string(), 
               parameters);
        A.set_frame(1); //Bjorken Frame

        /// Initialize a simple hadronizer
        JetDenseMediumHadronize Hadronizer;
        // Fill in all events
        std::vector<particle> plist, dplist, thermal_list, hlist;
        std::vector<current> clist;

        PGunWShower pythiagen(Q0, args["ic"].as<fs::path>().string());

        PGunWShower pythiagen2(0.5, args["ic"].as<fs::path>().string());
        JetFinder jf(300,300,3., need_response_table, args["response-table"].as<fs::path>().string());
        /// Initialzie a hydro reader
        std::vector<double> dE0({0,0,0,0,0}), dE({0,0,0,0,0});
        for (int i=0; i<args["nevents"].as<int>(); i++){
            plist.clear();
            clist.clear();
            pythiagen2.Generate(PID, args["E0"].as<double>(), plist);
            Hadronizer.hadronize(plist, hlist, thermal_list, Q0, Tf, 1);
            jf.MakeETower(0.6, Tf, 0.01, hlist, clist, 10, false);
            for (int ir=0; ir<5; ir++){
               
                 double rl = ir*.2, rh = (ir+1)*.2;
                 for (int ieta=0; ieta<300; ieta++){
                     double eta = -3+ieta*6./299.;
                     for (int iphi=0; iphi<300; iphi++){
                         double phi = -M_PI+iphi*2*M_PI/299.;
                         double R = std::sqrt(2.*(std::cosh(eta)
                                                 -std::cos(phi)));
                         if (rl <= R && R<rh) {
                             dE0[ir] += jf.PT[ieta][iphi];
                         }
                     }
                 }
            }
        }
          std::stringstream fheader1;
        fheader1 << args["output"].as<fs::path>().string() 
                 << processid << "-xJ.dat";
        std::ofstream f1(fheader1.str());

        for (int i=0; i<args["nevents"].as<int>(); i++){
            LOG_INFO << "___________________________";
            Medium<2> med1(args["hydro"].as<fs::path>().string());
            double mini_tau0 = med1.get_tauH();
            plist.clear();
            clist.clear();
            pythiagen.Generate(PID, args["E0"].as<double>(), plist);


             double Einitial = 0.;
             for (auto & p : plist){ 
                 auto p2 = p;
                 p2.p.boost_back(0., 0., std::tanh(p2.x.x3()));
                 double R = std::sqrt(2.*(std::cosh(p2.p.pseudorap())
                                          -std::cos(p2.p.phi())));
                 if (R<0.3) {
                     Einitial += p2.p.xT();
                 }
             }

            // move partciles below tau0 to tau0
            for (auto & p : plist){
                if (p.x.x0() >= mini_tau0 ) continue;
                else {
                    p.x.a[0] = mini_tau0;
                    p.x.a[1] += p.p.x()/p.p.t()*mini_tau0;
                    p.x.a[2] += p.p.y()/p.p.t()*mini_tau0;
                }
            }   
            std::vector<particle> new_plist, pOut_list;
            while(med1.load_next()) {
                new_plist.clear();
                double current_hydro_clock = med1.get_tauL();
                double dtau = med1.get_hydro_time_step();
                LOG_INFO << "Hydro t = " 
                         << current_hydro_clock/5.076 << " fm/c";
                for (auto & p : plist){     
                    if ( std::abs(p.x.x3())>5. ){
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
                        A.update_single_particle(DeltaTau, 
                                                 T, {vx, vy, vz}, 
                                                 p, pOut_list
                                                 );  
                        for (auto & fp : pOut_list) {
                            // compute energy momentum loss of hard partons
                            ploss = ploss - fp.p;
                            new_plist.push_back(fp);
                        }
                            current J; 
                            J.p = ploss;
                            J.etas = p.x.x3();
                            clist.push_back(J); 
                    }
                }
                plist = new_plist;
            }
            // back to lab frame
            for (auto & p : plist) {
                p.p.boost_back(0., 0., std::tanh(p.x.x3()));
                p.radlist.clear();
            }
            Hadronizer.hadronize(plist, hlist, thermal_list, Q0, Tf, 1);

                for(auto & it : thermal_list){
                    double vz = std::tanh(it.x.x3());
                    current J; 
                    J.p = it.p.boost_to(0, 0, vz)*(-1.);
                    J.etas = it.x.x3();
                    clist.push_back(J);  
                }
            jf.MakeETower(0.6, Tf, 0.01, hlist, clist, 10, false);

             double Efinal = 0.;
             for (int ieta=0; ieta<300; ieta++){
                     double eta = -3+ieta*6./299.;
                     for (int iphi=0; iphi<300; iphi++){
                         double phi = -M_PI+iphi*2*M_PI/299.;
                         double R = std::sqrt(2.*(std::cosh(eta)
                                                 -std::cos(phi)));
                         if (R<0.3) {
                            Efinal += jf.PT[ieta][iphi];
                         }
                     }
                 }
            f1 << Einitial << " " << Efinal << std::endl;
            
            for (int ir=0; ir<5; ir++){
               
                 double rl = ir*.2, rh = (ir+1)*.2;
                 for (int ieta=0; ieta<300; ieta++){
                     double eta = -3+ieta*6./299.;
                     for (int iphi=0; iphi<300; iphi++){
                         double phi = -M_PI+iphi*2*M_PI/299.;
                         double R = std::sqrt(2.*(std::cosh(eta)
                                                 -std::cos(phi)));
                         if (rl <= R && R<rh) {
                             dE[ir] += jf.PT[ieta][iphi];
                         }
                     }
                 }

            }
        }
        std::stringstream fheader2;
        fheader2 << args["output"].as<fs::path>().string() 
                 << processid << "-PT.dat";
        std::ofstream f2(fheader2.str());
        for (int ir=0; ir<5; ir++) 
            f2 << dE0[ir]/args["nevents"].as<int>() << " ";
        f2 << std::endl;
        for (int ir=0; ir<5; ir++) 
            f2 << dE[ir]/args["nevents"].as<int>() << " ";
        f2 << std::endl;
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



