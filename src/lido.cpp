#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include "Langevin.h"
#include "lido.h"
#include "random.h"
#include "simpleLogger.h"

lido::lido(std::string setting_path, std::string table_path, 
         std::vector<double> parameters){
    Lido_Ecut = parameters[3];
    std::string table_mode = (boost::filesystem::exists(table_path)) ? 
                         "old" : "new";
    CM = std::unique_ptr<collision_manager>(
           new collision_manager(table_mode, setting_path, table_path, parameters)
    );
}


fourvec ParallelTransportAlongZ(
          fourvec p, coordinate x, double dx3, int Frame){
    if (Frame==0) return p;
    else if (Frame==1){
        double ch = std::cosh(dx3), sh = std::sinh(dx3);
        double newp0 =  p.t()*ch - p.z()*sh; 
        double newpz = -p.t()*sh + p.z()*ch;
        return fourvec{newp0, p.x(), p.y(), newpz};
    }
}


void lido::FreeStream(
    particle & p, double dt){
    if (FrameChoice==0){ // Cartesian
        // p: E, px, py, pz in lab frame
        // x: t, x, y, z in lab frame
        double dtoverE = dt/p.p.t();
        p.x.a[0] = p.x.x0() + dt;
        p.x.a[1] = p.x.x1() + p.p.x()*dtoverE;
        p.x.a[2] = p.x.x2() + p.p.y()*dtoverE;
        p.x.a[3] = p.x.x3() + p.p.z()*dtoverE;
    }
    else if (FrameChoice==1){ // Bjorken
        // dt is dtau
        double dtau = dt;
        // x is tau, x, y, etas
        // p is the four moemtnum (E', px', py', pz' [GeV]) in frame vF = (0, 0, z/t=tanh(etas))
        // Apply inertia force for dtau/2
        double m2 = p.mass*p.mass;
        p.p.a[3] -= dtau/2.*p.p.a[3]/p.x.x0();
        p.p.a[0] = std::sqrt(m2 + p.p.pabs2());
        // Update x to t+dtau
        double dt_over_Eprime = dtau/p.p.t();
        p.x.a[0] += dtau;
        p.x.a[1] += p.p.x()*dt_over_Eprime;
        p.x.a[2] += p.p.y()*dt_over_Eprime;
        p.x.a[3] += p.p.z()/p.x.x0()*dt_over_Eprime;
        // Apply inertia force for dtau/2
        p.p.a[3] -= dtau/2.*p.p.a[3]/p.x.x0();
        p.p.a[0] = std::sqrt(m2 + p.p.pabs2());
    }
}

void lido::Diffusion(
    particle & p, double dt, double T, std::vector<double> v){
    // In the frame associated to the coordiantes
    // e.g. in Bjorken frame dt=dtau, and v is the medium velocity
    //                  in the frame vF = (0, 0, z/t=tanh(etas));
    auto pcell = p.p.boost_to(v[0], v[1], v[2]);
    double dtcell = dt*pcell.t()/p.p.t();
    fourvec pnew;
    for (int i=0; i<5; i++){
        Ito_update_rest(p.pid, dtcell/5., p.mass, T, pcell, pnew);
        if (p.is_virtual) pcell = pnew*(pcell.t()/pnew.t());
        else pcell = pnew;
        pcell = put_on_shell(pcell, p.pid);
    }
    p.p = pcell.boost_back(v[0], v[1], v[2]);
}

int lido::update_single_particle(
    double dt, double T, std::vector<double> v3,
    particle & pIn, std::vector<particle> & pOut_list){
    // In the frame associated to the coordiantes
    // e.g. in Bjorken frame dt=dtau, and v is the medium velocity
    //                  in the frame vF = (0, 0, z/t=tanh(etas));
    pOut_list.clear();
    int abspid = std::abs(pIn.pid);
    pIn.vcell[0] = v3[0];
    pIn.vcell[1] = v3[1];
    pIn.vcell[2] = v3[2];
    pIn.Tf = T;
    //---A---: FreeStream the particle
    FreeStream(pIn, dt);


    // we only handle u,d,s,c,b,g
    // other wise, leave it untouched
    bool isParton = (abspid==21) || (abspid<=5);
    bool isLightParton = (abspid==21) || (abspid<=3);
    bool isHeavyQuark = (abspid==4) || (abspid==5);
    if ( !isParton ){
        pIn.radlist.clear();
        pOut_list.push_back(pIn);
        return pOut_list.size();
    }

    double Emin = Lido_Ecut * T;
    double mD2 = t_channel_mD2->get_mD2(T);
    double mD = std::sqrt(mD2);
    double minf = mD/std::sqrt(2);
    double Eradmin = 3.*T;

    // If a light parton has
    // E<Emin and Q<mD
    // and is not a virtual particle, drop it from LIDO
    double Ecell = pIn.p.boost_to(v3[0], v3[1], v3[2]).t();
    if ( (!pIn.is_virtual) && (pIn.Q0<minf) && (Ecell < Emin) && isLightParton
    ){
        pIn.radlist.clear();
        return pOut_list.size();
    }
 
    //---B---: If the pareton is not yet formed,
    // skip its interactions
    bool formed_from_vac = (pIn.x.x0()>=pIn.tau0);
    if (!formed_from_vac) {
        pIn.radlist.clear();
        pOut_list.push_back(pIn);
        return pOut_list.size();
    }

    //---C---: Diffusion
    Diffusion(pIn, dt, T, v3);

    //---D---: Sample scatterings
    std::vector<particle> FS;
    int process_id = CM->sample(pIn, FS, T, v3, dt);

    //---E---: return virtual particle copy
    if (pIn.is_virtual){  
        // virtual particle only evolves according to semi-classical
        // 2-2 collisions, only return itself
        if (process_id>0) pIn.p = FS[0].p;
        pIn.radlist.clear();
        pOut_list.push_back(pIn);
        return pOut_list.size();
    }
    
    //---F---: Handle scattering final states
    // From now on, pIn is a real particle
    if (process_id==-1){
        // Nothing happens      
    }
    else if (process_id==-2){
        // Something happened, but sampling failed, 
        // issue a warining, continue
        LOG_WARNING << "Sampling failed, return the initial hard parton";
    }
    else if (is2to2(process_id)){
        // Two-body collision
        // add both final-state particle to Lido list
        assign_2to2_color(process_id, 
                       pIn.pid, FS[0].pid, FS[1].pid,
                       pIn.col, pIn.acol,
                       FS[0].col, FS[0].acol,
                       FS[1].col, FS[1].acol
                       );
        if (isPairProduction(process_id)) {
            // 2->2 pair production destroy all the interference
            pIn.radlist.clear();
            pOut_list.push_back(FS[0]);
            pOut_list.push_back(FS[1]);
            return pOut_list.size();
        }
        else {
            // semi-classical 2->2 preserve the identity
            pIn.p = FS[0].p;
            pIn.col = FS[0].col;
            pIn.acol = FS[0].acol;
            // add the recoil particle
            pOut_list.push_back(FS[1]);
        }
    }
    else if (is2to2_nlo(process_id)){
        // elastic for hard 
        // soft parton undergoes quark-pair production
        assign_2to2_nlo_color(process_id, 
                       pIn.pid, FS[0].pid, FS[1].pid, FS[2].pid,
                       pIn.col, pIn.acol,
                       FS[0].col, FS[0].acol,
                       FS[1].col, FS[1].acol,
                       FS[2].col, FS[2].acol
                       );
        // hard parton changes momentum and color
        pIn.p = FS[0].p;
        pIn.col = FS[0].col;
        pIn.acol = FS[0].acol;
        // add two recoil particle
        pOut_list.push_back(FS[1]);
        pOut_list.push_back(FS[2]);
    }
    else if (is1to2(process_id)){
        // if the radiated parton is energytic enough, take it as a 
        // virtual particle about to from
        if (  (FS[0].p.boost_to(v3[0], v3[1], v3[2]).t() > Eradmin || std::abs(FS[0].pid)==4 || std::abs(FS[0].pid)==5 )
           && (FS[1].p.boost_to(v3[0], v3[1], v3[2]).t() > Eradmin || std::abs(FS[1].pid)==4 || std::abs(FS[1].pid)==5 )
        ){
            FS[0].is_virtual = true;
            FS[0].T0 = T;
            FS[1].is_virtual = true;
            FS[1].T0 = T;
            pIn.radlist.push_back({FS[0], FS[1]});
        }
    }
    else if (is2to3(process_id)){
        // if the radiated parton is energytic enough, take it as a 
        // virtual particle about to from
        if (  (FS[0].p.boost_to(v3[0], v3[1], v3[2]).t() > Eradmin || std::abs(FS[0].pid)==4 || std::abs(FS[0].pid)==5 )
           && (FS[2].p.boost_to(v3[0], v3[1], v3[2]).t() > Eradmin || std::abs(FS[2].pid)==4 || std::abs(FS[2].pid)==5 )
        ){
            FS[0].is_virtual = true;
            FS[0].T0 =T;
            FS[2].is_virtual = true;
            FS[0].T0 =T;
            pIn.radlist.push_back({FS[0], FS[2]});
        }
        // recoil parton is neglected
    }

    //---G---: LPM correction
    // Only for real particle that has a non-empty radiation list

    if ((!pIn.radlist.empty()) && (!pIn.is_virtual)){
        // Remove those rad particle whose kinematic constrain is  
        // violated
        for (std::vector<std::vector<particle> >::iterator 
             it = pIn.radlist.begin(); 
             it != pIn.radlist.end();){
            fourvec PB = ParallelTransportAlongZ(
                (*it)[1].p, (*it)[1].x0, 
                pIn.x.x3()-(*it)[1].x.x3(), FrameChoice
            );
            fourvec PC = ParallelTransportAlongZ(
                (*it)[0].p, (*it)[0].x0, 
                pIn.x.x3()-(*it)[0].x.x3(), FrameChoice
            );
            if(PB.t() + PC.t() >= pIn.p.t() ) 
                it = pIn.radlist.erase(it);
            else it++;
        }

        // loop over each virtual particles and evolve it
        for (std::vector<std::vector<particle> >::iterator 
             it = pIn.radlist.begin(); 
             it != pIn.radlist.end(); it++){
            std::vector<particle> F1, F2;
            update_single_particle(dt, T, v3, (*it)[0], F1);
            (*it)[0].p = F1[0].p;
            update_single_particle(dt, T, v3, (*it)[1], F2);
            (*it)[1].p = F2[0].p;
        }
        bool nonclassical = false;
        // check which one has reached its formation time
        for (std::vector<std::vector<particle> >::iterator 
             it = pIn.radlist.begin(); 
             it != pIn.radlist.end();){
            if (nonclassical) break;
            double taun, kt2n;
            // to compute formation time in Bjorken frame,
            // we will need to parallel transport other vectors 
            // to the coordinate of the current partcile
            // 1) Daughter partons: boost in eta
            // 2) Original copy: transport in tau
            fourvec P0 = (*it)[1].mother_p.boost_back(
                           0, 0, std::tanh((*it)[1].x0.x3())).boost_to(
                           0, 0, std::tanh(pIn.x.x3()));
            fourvec PB = ParallelTransportAlongZ(
                (*it)[1].p, (*it)[1].x0, 
                pIn.x.x3()-(*it)[1].x.x3(), FrameChoice
            );
            fourvec PB0 = ParallelTransportAlongZ(
                (*it)[1].p0, (*it)[1].x0, 
                pIn.x.x3()-(*it)[1].x.x3(), FrameChoice
            );
            fourvec PC = ParallelTransportAlongZ(
                (*it)[0].p, (*it)[0].x0, 
                pIn.x.x3()-(*it)[0].x.x3(), FrameChoice
            );
           
            formation_time(taun, kt2n, 
                          pIn.pid, (*it)[1].pid, (*it)[0].pid, 
                          pIn.p, PB, PC, T, P0);
            // it it has not formed, continue
            if (pIn.x.x0()-(*it)[1].x0.x0() <= taun) it++;
            else{
                double P0abs = P0.pabs();
                fourvec nbar{1., -P0.x()/P0abs, -P0.y()/P0abs, -P0.z()/P0abs};
                double Eplus_mother = dot(P0, nbar);
                double Z = dot(PB, nbar)/Eplus_mother;
                // if formed, apply LPM suppression in the Bjorken frame
                double Acceptance = 0.;
                double kt20 = measure_perp(P0, PB0).pabs2();
                double kt21 = measure_perp(P0, PB).pabs2();
                double kt22 = measure_perp(P0, PC).pabs2();
                double Running = std::min(
                       alpha_s(kt21,T)/alpha_s(kt20,T), 1.);
                
                // an NLL-inspired suppression factor for the LPM effect
                double mD2 = t_channel_mD2->get_mD2(T);
                double lnQ2_1 = std::log(1. + taun/(*it)[1].mfp0);
                double lnQ2_0 = std::log(1. + 
                 6.*Z*(1-Z)*pIn.p.boost_to(v3[0], v3[1], v3[2]).t()
                   * (*it)[1].T0 / mD2);
                double MassFactors = 1.0;
                if (std::abs(pIn.pid)==4 || std::abs(pIn.pid)==5){
                    double m2 = std::pow(Z*pIn.mass, 2);
                    MassFactors = std::pow(kt21/(kt21+m2), 2);
                }
                if (std::abs(pIn.pid)==21 && (*it)[1].pid!=21){
                    double m2 = std::pow((*it)[1].mass, 2);
                    MassFactors = std::pow(kt21/(kt21+m2), 2);
                }
                double log_factor = std::sqrt(lnQ2_1/lnQ2_0);
                double LPM = std::min((*it)[1].mfp0 / taun * log_factor, 1.);
                // The final acceptance factor
                Acceptance = LPM * Running * MassFactors;

                if (Srandom::rejection(Srandom::gen) < Acceptance){   
                    assign_n2np1_color(
                       pIn.pid, (*it)[0].pid, (*it)[1].pid,
                       pIn.col, pIn.acol,
                       (*it)[0].col, (*it)[0].acol,
                       (*it)[1].col, (*it)[1].acol
                       );
                    if ((*it)[1].pid==21) { 
                        pIn.p = pIn.p - PB;
                    }
                    if ((*it)[1].pid!=21) {
                        pIn.p = pIn.p - PB;   
                        pIn.pid = -(*it)[1].pid;
                        pIn.mass = pid2mass(pIn.pid);
                        nonclassical = true;
                    }
                    pIn.p = put_on_shell(pIn.p, pIn.pid);;
                    pIn.p0 = pIn.p; 

                    (*it)[1].p = put_on_shell(PB, (*it)[1].pid);
                    (*it)[1].p0 = (*it)[1].p;  
                    (*it)[1].x = pIn.x;
                    (*it)[1].x0 = pIn.x;

                    double Scale = std::sqrt(kt2n+mD2);
                    pIn.Q0 = Scale;
                    pIn.Q00 = Scale;
                 
                    pIn.col = (*it)[0].col;
                    pIn.acol = (*it)[0].acol;
                    (*it)[1].Q0 = Scale;
                    (*it)[1].Q00 = Scale;

                    (*it)[1].is_virtual = false;
                    pOut_list.push_back((*it)[1]);
                }
                // final remove it from the radlist
                it = pIn.radlist.erase(it);
            }
        }
        if (nonclassical) pIn.radlist.clear();
    }



    // Add the mother parton to the end of the output list
    //also transform back its pid
    pOut_list.push_back(pIn);
    return pOut_list.size();
}


