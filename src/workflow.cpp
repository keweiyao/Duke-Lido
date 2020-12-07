#include "workflow.h"
#include "simpleLogger.h"
#include <fstream>
#include "matrix_elements.h"
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/variant/get.hpp>
#include <boost/any.hpp>
#include <boost/foreach.hpp>
#include "logo.h"
#include "Langevin.h"
#include "predefine.h"
#include "random.h"



void initialize(std::string mode, std::string setting_path, std::string table_path,
		double mu, double afix, double theta, double cut){
    print_logo();
    boost::property_tree::ptree config;
    std::ifstream input(setting_path);
    read_xml(input, config);
    // 0: constant time step, 1: constant proper time step
    time_type = config.get("t_type", 0); 
    // 0: default, 1: local approxiamtion
    bool Adiabatic_LPM = config.get("LPM_type", 0);
}

int update_particle_momentum_Lido(
    double dt_input, double temp, std::vector<double> v3cell,
    particle & pIn, std::vector<particle> & pOut_list, double Tf){
    double Emin = Lido_Ecut * temp;
        // minimum rad energy
    auto p00 = pIn.p;
    pOut_list.clear();
    pIn.Tf = temp;
    pIn.vcell[0] = v3cell[0];
    pIn.vcell[1] = v3cell[1];
    pIn.vcell[2] = v3cell[2];
    auto x0 = pIn.x;
    double dt_for_pIn = compute_realtime_to_propagate(dt_input, pIn.x, pIn.p);
    pIn.freestream(dt_for_pIn);
    // we only handle u,d,s,c,b,g
    bool formed_from_vac = (pIn.x.t()>=pIn.tau_i);
    if ( (!((pIn.pid=21) || (std::abs(pIn.pid)<=5)) )
       || (temp<Tf) || (!formed_from_vac) ){
        pIn.radlist.clear();
        pOut_list.push_back(pIn);
        return pOut_list.size();
    }

    double mD2 = t_channel_mD2->get_mD2(temp);
    double mD = std::sqrt(mD2);
    double minf = mD/1.414;
    double Eradmin = minf;

    // Don't touch particles already below soft cut
    double E0_cell = pIn.p.boost_to(v3cell[0], v3cell[1], v3cell[2]).t();
    if ((!pIn.is_virtual) && (pIn.Q0<minf)
        && (E0_cell < Emin)
        && (std::abs(pIn.pid)==1 || std::abs(pIn.pid)==2 || 
            std::abs(pIn.pid)==3 || std::abs(pIn.pid)==21) 
	){
        pIn.radlist.clear();
        return pOut_list.size();
    }


    int channel_pid;
    if (std::abs(pIn.pid) <= 3) channel_pid = 123;
    else channel_pid = std::abs(pIn.pid);

    if (E0_cell<Emin){
        pIn.radlist.clear();
        pOut_list.push_back(pIn);
        return pOut_list.size();
    }


    for (int i=0; i<10; i++){
        fourvec pnew;
        Ito_update(pIn.pid, dt_for_pIn/10, pIn.mass, temp, v3cell, pIn.p, pnew, pIn.is_virtual);
	if (pIn.is_virtual) pIn.p = pnew*(pIn.p.t()/pnew.t());
	else pIn.p = pnew;
        pIn.p.a[0] = std::sqrt(pIn.p.pabs2()+pIn.mass*pIn.mass);
    }



    // If a scattering happens:
    if (channel >= 0){


        // elastic process changes the momentum immediately
        if (channel == 0 || channel == 1){
            pIn.p = FS[0];
            if (!pIn.is_virtual){  
		    int id = (channel==0) ? (Srandom::sample_flavor(3)) : 21;
		    int col=-100, acol=-100, mcol=-100, macol=-100;
		    SampleFlavorAndColor(
		        pIn.pid, pIn.col, pIn.acol, channel,
		        id, col, acol, mcol, macol);
		    particle ep = produce_parton(id, FS[1], pIn, false);
		    ep.origin = 2;
		    ep.col = col; 
		    ep.acol = acol;
		    pIn.col = mcol;
		    pIn.acol = macol;
		    pOut_list.push_back(ep);
           }
        }
        // inelastic process takes a finite time to happen
        if (channel == 2 || channel == 3){
            // add the radiating daughters as ``virtual" particles
            // these are virtual particles, that pops out at a rate
            // given by the incoherent calculation, which will be
            // dropped (suppressed) according to the LPM effect
            // only handle > mD particles
            if (FS.size()==3){
            if (FS[2].boost_to(v3cell[0], v3cell[1], v3cell[2]).t() > Eradmin) { 
                particle vp = produce_parton(21, FS[2], pIn, true);
                // The local 2->2 mean-free-path is estimated with
                // the qhat_hard integrate from the 2->2 rate
                double xfrac = vp.p.t() / pIn.p.t();
                double local_qhat = 0.;
                double vp_cell = vp.p.boost_to(v3cell[0], v3cell[1], v3cell[2]).t();
                double Lnvp_cell = std::log(vp_cell);
                double boost_factor = vp.p.t() / vp_cell;
                BOOST_FOREACH (Process &r, AllProcesses[21]){
                    switch (r.which()){
                    case 0:
                        if (boost::get<Rate22>(r).IsActive()){
                            tensor A = boost::get<Rate22>(r).GetSecondM({Lnvp_cell, temp});
                            local_qhat += A.T[1][1] + A.T[2][2];
                        }
                        break;
                    default:
                        break;
                    }
                }  
                // estimate mfp in the lab frame
                vp.mfp0 = LPM_prefactor * mD2 / local_qhat * boost_factor;
                pIn.radlist.push_back(vp);
            }
            }
        }
        if (channel == 4 || channel == 5){
            // Absorption processes happens mostly for gluon energy ~ 3*T, therefore we negelected the LPM effect
            pIn.p = FS[0];
        }
        if (channel == 6){
            if (FS.size()==2){
            if (FS[1].boost_to(v3cell[0], v3cell[1], v3cell[2]).t() > Eradmin) {
                particle vp = produce_parton(21, FS[1], pIn, true);
                // estimate mfp in the lab frame
                vp.mfp0 = LPM_prefactor * mD2 / qhatg * pIn.p.t() / E_cell;
                pIn.radlist.push_back(vp);
            }
            }
        }
        if (channel == 7){
            // Absorption processes happens mostly for gluon energy ~ 3*T, therefore we negelected the LPM effect
            pIn.p = FS[0];
        }
        // elastic indued charm pair production
        if (channel == 8){
            // Two hard particle in final state;
            // Change of pid! clear all radiaiton 
            pIn.radlist.clear();
            pIn.p = FS[0];
            pIn.pid = 4;
	    int id = -4;
            int col=-100, acol=-100, mcol=-100, macol=-100;
	    SampleFlavorAndColor(
		pIn.pid, pIn.col, pIn.acol, channel,
		id, col, acol, mcol, macol);
	    particle ep = produce_parton(id, FS[1], pIn, false);
	    ep.origin = 2;
	    ep.col = col; 
	    ep.acol = acol;
	    pIn.col = mcol;
	    pIn.acol = macol;
            pOut_list.push_back(ep);
        }
        if (channel == 9 || channel == 10){
            if (FS.size()==3 && FS[0].t()>1.4 && FS[2].t()>1.4){
	        particle vp = produce_parton(Srandom::binary_choice()? 4:-4, FS[2], pIn, true);
                vp.mass = 1.3;
                vp.p = vp.p*(std::sqrt(vp.p.t()*vp.p.t()-vp.mass*vp.mass)
                           /vp.p.t());
                vp.p.a[0] = std::sqrt(vp.mass*vp.mass+vp.p.pabs2());
	        // The local 2->2 mean-free-path is estimated with
	        // the qhat_hard integrate from the 2->2 rate
	        double xfrac = vp.p.t() / pIn.p.t();
	        double local_qhat = 0.;
	        double vp_cell = vp.p.boost_to(v3cell[0], v3cell[1], v3cell[2]).t();
	        double Lnvp_cell = std::log(vp_cell);
	        double boost_factor = vp.p.t() / vp_cell;
	        BOOST_FOREACH (Process &r, AllProcesses[21]){
	            switch (r.which()){
	            case 0:
	                if (boost::get<Rate22>(r).IsActive()){
	                    tensor A = boost::get<Rate22>(r).GetSecondM({Lnvp_cell, temp});
	                    local_qhat += A.T[1][1] + A.T[2][2];
	                }
	                break;
	            default:
	                break;
	            }
	        }  
	        // estimate mfp in the lab frame
	        vp.mfp0 = LPM_prefactor * mD2 / local_qhat * boost_factor;
	        pIn.radlist.push_back(vp);
            }
        }
        if (channel == 11 && FS[0].t()>1.4 && FS[1].t()>1.4){
            if (FS.size()==2){
	        particle vp = produce_parton(
                          Srandom::binary_choice()? 4:-4, FS[1], pIn, true);
                vp.mass = 1.3;
                vp.p = vp.p*(std::sqrt(vp.p.t()*vp.p.t()-vp.mass*vp.mass)
                           /vp.p.t());
                vp.p.a[0] = std::sqrt(vp.mass*vp.mass+vp.p.pabs2());
                // estimate mfp in the lab frame
                vp.mfp0 = LPM_prefactor * mD2 / qhatg * pIn.p.t() / E_cell;
                pIn.radlist.push_back(vp);
            }
        }
    }

    bool has_altered_ID = false;
    // Handle the virtual particle, (only for real mother parton)
    if ((!pIn.radlist.empty()) && (!pIn.is_virtual)){
        // loop over each virtual particles
        for (std::vector<particle>::iterator it = pIn.radlist.begin(); it != pIn.radlist.end();){
            int split_type = 0;
            int abs_rad_id = std::abs(it->pid);
            if (pIn.pid != 21 &&  abs_rad_id==21) split_type = 1; // q --> q + g
            if (pIn.pid == 21){
                if (abs_rad_id == 21) split_type = 2; // g --> g + g
                else if (abs_rad_id<4) split_type = 3; // g --> q + qbar
                else if (abs_rad_id==4 || abs_rad_id==5) split_type = 4; // g --> Q + Qbar
                else LOG_FATAL << "Wrong splitting!";
            }
            double taun = formation_time(it->mother_p, it->p, temp, split_type);
            std::vector<particle> pnew_Out;

            // In the local (adiabatic LPM) mode,
            // use local information to determine the number of rescatterings
            if (Adiabatic_LPM){
                do{
                   update_particle_momentum_Lido(dt_input, temp, v3cell, *it, pnew_Out, Tf);
                   taun = formation_time(it->mother_p, it->p, temp, split_type);
                }while(it->x.t() - it->x0.t() <= taun);
                // once finished, reset its x to where it is produced
                it->x = it->x0;
            }
            // update the partilce for this time step
            update_particle_momentum_Lido(
                           dt_input, temp, v3cell, *it, pnew_Out, Tf);
            if (it->x.t() - it->x0.t() <= taun && !Adiabatic_LPM ) {
                it++;
            }
            else{
                // if formed, apply LPM suppression
                double Acceptance = 0.;
                // for medium-induced radiation
                // 1): a change of running-coupling from elastic broadening
                double kt20 = measure_perp(it->mother_p, it->p0).pabs2();
                double kt2n = measure_perp(pIn.p, it->p).pabs2();
                double Running = std::min(
			alpha_s(kt2n,it->T0)/alpha_s(kt20,it->T0), 1.);
                // 2): a dead-cone approximation for massive particles

		double DeadCone = 1.;
                if (split_type==1){
                    double theta2 = kt2n / std::pow(it->p.t(), 2);
                    double thetaM2 = std::pow(pIn.mass / it->mother_p.t(), 2);
                    DeadCone = std::pow(theta2/(theta2 + thetaM2), 2);
                }
                if (split_type==4) 
                    DeadCone = std::pow(kt2n/(kt2n+std::pow(it->mass,2)), 2);
         
                // 3): an NLL-inspired suppression factor for the LPM effect
                double mD2 = t_channel_mD2->get_mD2(temp);
                double lnQ2_1 = std::log(1. + taun / it->mfp0);
                double lnQ2_0 = std::log(1. + 6 * it->p.boost_to(v3cell[0], v3cell[1], v3cell[2]).t() * temp / mD2);
                double log_factor = std::sqrt(lnQ2_1 / lnQ2_0);
                double LPM = std::min(it->mfp0 / taun * log_factor, 1.);
                // The final rejection factor
                Acceptance = LPM * Running * DeadCone;

                if (Srandom::rejection(Srandom::gen) < Acceptance){
                    // accepted branching causes physical effects
                    // momentum change, and put back on shell
                    if (split_type==1 || split_type==2) {
                        // mother parton id unchanged
                        pIn.p = pIn.p - it->p;
                        it->p0 = it->p;                 
                        pIn.p.a[0] = std::sqrt(pIn.mass*pIn.mass+pIn.p.pabs2());
		        pIn.p0 = pIn.p;
                        int col=-100, acol=-100, mcol=-100, macol=-100;
                        SampleFlavorAndColor(pIn.pid, pIn.col, pIn.acol, 2, it->pid, 
                                            col, acol, mcol, macol);
                        it->col = col; 
                        it->acol = acol;
                        pIn.col = mcol;
                        pIn.acol = macol;
                    }
                    if (split_type==3 || split_type==4) {
                        // mother parton id unchanged
                        it->mass = 1.3;
                        pIn.mass = it->mass;
                        pIn.pid = -it->pid;
                        pIn.p = pIn.p - it->p;
                        it->p0 = it->p;                 
                        pIn.p.a[0] = std::sqrt(pIn.mass*pIn.mass+pIn.p.pabs2());
                        it->p.a[0] = std::sqrt(it->mass*it->mass+it->p.pabs2());
		        pIn.p0 = pIn.p;
                        int col=-100, acol=-100, mcol=-100, macol=-100;
                        SampleFlavorAndColor(pIn.pid, pIn.col, pIn.acol, 
                                            9, it->pid, 
                                            col, acol, mcol, macol);
                        it->col = col; 
                        it->acol = acol;
                        pIn.col = mcol;
                        pIn.acol = macol;
                        has_altered_ID = true;
                    }
                    double Scale = std::sqrt(kt2n+mD2);
		    pIn.Q0 = Scale;
                    pIn.Q00 = Scale;
	            it->Q0 = Scale;
		    it->Q00 = Scale;
                    it->is_virtual = false;
		    if (pIn.origin == 2) it->origin=2;
                    if (pIn.origin == 0 || pIn.origin==1) it->origin=1;
                    pOut_list.push_back(*it);
                }
                // remove it from the radlist
                it = pIn.radlist.erase(it);
            }
        }
    }
    if (has_altered_ID) pIn.radlist.clear();

    for (std::vector<particle>::iterator it = pIn.radlist.begin(); it != pIn.radlist.end();){
            if(it->p.t()>=pIn.p.t()) it = pIn.radlist.erase(it);
            else it++;
    }


    // Add the mother parton to the end of the output list
    //also transform back its pid
    pOut_list.push_back(pIn);

    // Compute energy source term using In and Out states of the hard-parton system
    return pOut_list.size();
}


void output_oscar(const std::vector<particle> plist,
                  int abspid, std::string fname){
    // output OSCAR Format
    int Nparticles = 0;
    for (auto &p : plist)
    {
        if (std::abs(p.pid) == abspid)
            Nparticles++;
    }

    std::ofstream f(fname);
    f << "OSC1997A\n";
    f << "final_id_p_x\n";
    f << "    lbt  1.0alpha   208    82   208    82   aacm  0.1380E+04        1\n";
    f << "        1  " << std::setw(10)
      << Nparticles << "    0.001    0.001    1    1       1\n";

    int i = 0;
    for (auto &p : plist){
        if (std::abs(p.pid) == abspid){
            f << std::setw(10) << std::setfill(' ') << i << "  "      // particle index, i10,2x
              << std::setw(10) << std::setfill(' ') << p.pid << "  "; // particle id, i10,2x
            f << ff(p.p.x()) << "  "
              << ff(p.p.y()) << "  "
              << ff(p.p.z()) << "  "
              << ff(p.p.t()) << "  "
              << ff(p.mass) << "  "
              << ff(p.x.x() / fmc_to_GeV_m1) << "  "
              << ff(p.x.y() / fmc_to_GeV_m1) << "  "
              << ff(p.x.z() / fmc_to_GeV_m1) << "  "
              << ff(p.x.t() / fmc_to_GeV_m1) << "  "
              << ff(p.Tf) << "  "
              << ff(p.vcell[0]) << "  "
              << ff(p.vcell[1]) << "  "
              << ff(p.vcell[2]) << "  "
              << ff(p.p0.x()) << "  "
              << ff(p.p0.y()) << "  "
              << ff(p.p0.z()) << "  "
              << ff(p.weight) << "\n";
        }
        i++;
    }
}
