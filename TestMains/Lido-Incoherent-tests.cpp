#include <string>
#include <iostream>
#include <exception>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/program_options.hpp>

#include "simpleLogger.h"
#include "lido.h"
namespace po = boost::program_options;
namespace fs = boost::filesystem;

int main(int argc, char* argv[]){
    using OptDesc = po::options_description;
    OptDesc options{};
    options.add_options()
          ("help,h", "show this help message and exit")
          ("lido-setting,s", 
            po::value<fs::path>()->value_name("PATH")->required(),
           "Lido table setting file")
          ("lido-table,t", 
            po::value<fs::path>()->value_name("PATH")->required(),
           "Lido table path to file")
	  ("muT",
           po::value<double>()->value_name("DOUBLE")->default_value(1.5,"1.5"),
	   "mu_min/piT")
	  ("afix",
           po::value<double>()->value_name("DOUBLE")->default_value(-1.,"-1."),
           "fixed alpha_s, <0 for running alphas")
           ("cut",
           po::value<double>()->value_name("DOUBLE")->default_value(4.,"4."),
           "cut between diffusion and scattering, Qc^2 = cut*mD^2")
         ;
    po::variables_map args{};
    try{
        po::store(po::command_line_parser(argc, argv).options(options).run(), args);
        if (args.count("help")){
                std::cout << "usage: " << argv[0] << " [options]\n"
                          << options;
                return 0;
        }    
        if (!args.count("lido-setting")){
            throw po::required_option{"<lido-setting>"};
            return 1;
        }
        else {
            if (!fs::exists(args["lido-setting"].as<fs::path>())){
                throw po::error{"<lido-setting> path does not exist"};
                return 1;
            }
        }
        if (!args.count("lido-table")){
            throw po::required_option{"<lido-table>"};
            return 1;
        }
        double muT = args["muT"].as<double>();
        double cut = args["cut"].as<double>();
        double afix = args["afix"].as<double>();
        std::vector<double> parameters{muT, afix, cut, 4.};
        lido A(args["lido-setting"].as<fs::path>().string(), 
            args["lido-table"].as<fs::path>().string(), 
            parameters);
        A.set_frame(0);

        std::vector<particle> plist, pOut_list, newplist;
        double T0 = 0.4, dt = .02*5.076;
        int N = 400000; 
        fourvec p0{100,100,0,0};
        coordinate x0{0,0,0,0};
        plist.clear();
        int pid0 = 21;
        for (int i=0; i<N; i++){
            plist.push_back(make_parton(pid0, 0, 0, 0, p0, x0) );
        }
        for (int i=0; i<500; i++){
            double T = T0;//*std::pow(.5/(dt*i/5.076+.5), 1./3.);
            LOG_INFO << i << "N = " << plist.size();
            for (auto & p : plist){
                A.update_single_particle(dt, T, {0,0,0}, p, pOut_list);
               
                //for (auto & ip: pOut_list) newplist.push_back(ip);
                p.p = p.p*(p0.t()/p.p.t());
                p.pid = pid0;
                p.p = put_on_shell(p.p, p.pid);
            }
            double E0=0;
            for (auto &p:plist) E0+=p.p.t();
            LOG_INFO <<"dE/dL-----------: "<< (p0.t()-E0/N)/(i+1)/dt;
            //plist = newplist;
            //newplist.clear();
            
        }
        //for (auto & p : plist){
        //     f1 << p.pid << " " << p.p << " " << p.x << std::endl;
        //}

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


