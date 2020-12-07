#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include "Langevin.h"
#include "lido.h"
#include "random.h"
#include "simpleLogger.h"

lido::lido(std::string setting_path, std::string table_path, 
         std::vector<double> parameters){
    Lido_Ecut = parameters[2];
    std::string table_mode = (boost::filesystem::exists(table_path)) ? 
                         "old" : "new";
    CM = std::unique_ptr<collision_manager>(
           new collision_manager(table_mode, setting_path, table_path, parameters)
    );
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
    for (int i=0; i<10; i++){
        Ito_update_rest(p.pid, dtcell/10, p.mass, T, pcell, pnew);
        if (p.is_virtual) pcell = pnew*(pcell.t()/pnew.t());
        else pcell = pnew;
        pcell = put_on_shell(pcell, p.pid);
    }
    p.p = pcell.boost_back(v[0], v[1], v[2]);
}
std::ofstream f("stat.dat");
int lido::update_single_particle(
    double dt, double T, std::vector<double> v3,
    particle & pIn, std::vector<particle> & pOut_list){
    // In the frame associated to the coordiantes
    // e.g. in Bjorken frame dt=dtau, and v is the medium velocity
    //                  in the frame vF = (0, 0, z/t=tanh(etas));
    pOut_list.clear();
    int abspid = std::abs(pIn.pid);
    pIn.Tf = T;
    pIn.vcell[0] = v3[0];
    pIn.vcell[1] = v3[1];
    pIn.vcell[2] = v3[2];

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
    double Eradmin = minf;

    // If a light parton has
    // E<Emin and Q<mD
    // and is not a virtual particle, drop it from LIDO
    double Ecell = pIn.p.boost_to(v3[0], v3[1], v3[2]).t();
    if ( (!pIn.is_virtual) && (pIn.Q0<mD) && (Ecell < Emin) && isLightParton
    ){
        pIn.radlist.clear();
        return pOut_list.size();
    }
   

    //---A---: FreeStream the particle
    FreeStream(pIn, dt);
 
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
            // add the recoil particle
            pOut_list.push_back(FS[1]);
        }
    }
    else if (is1to2(process_id)){
        // if the radiated parton is energytic enough, take it as a 
        // virtual particle about to from
        if (FS[1].p.boost_to(v3[0], v3[1], v3[2]).t() > Eradmin){
            FS[1].is_virtual = true;
            pIn.radlist.push_back(FS[1]);
        }
    }
    else if (is2to3(process_id)){
        // if the radiated parton is energytic enough, take it as a 
        // virtual particle about to from
        if (FS[2].p.boost_to(v3[0], v3[1], v3[2]).t() > Eradmin){
            FS[2].is_virtual = true;
            pIn.radlist.push_back(FS[2]);
        }
        // recoil parton is neglected
    }

    //---G---: LPM correction
    // Only for real particle that has a non-empty radiation list
    
    if ((!pIn.radlist.empty()) && (!pIn.is_virtual)){
        // loop over each virtual particles and evolve it
        for (std::vector<particle>::iterator it = pIn.radlist.begin(); 
             it != pIn.radlist.end(); it++){
            std::vector<particle> virtual_FS;
            update_single_particle(dt, T, v3, *it, virtual_FS);
            it->p = virtual_FS[0].p;
        }
        bool nonclassical = false;
        // check which one has reached its formation time
        for (std::vector<particle>::iterator it = pIn.radlist.begin(); 
             it != pIn.radlist.end();){
            if (nonclassical) break;
            double taun = formation_time(it->mother_p, it->p, 
                                         T, pIn.pid, it->pid);
            // it it has not formed, continue
            if (it->x.x0()-it->x0.x0() <= taun) it++;
            else{
                // if formed, apply LPM suppression
                double Acceptance = 0.;
                double kt20 = measure_perp(it->mother_p, it->p0).pabs2();
                double kt2n = measure_perp(pIn.p, it->p).pabs2();
                double Running = std::min(
                       alpha_s(kt2n,T)/alpha_s(kt20,T), 1.);

                // an NLL-inspired suppression factor for the LPM effect
                double mD2 = t_channel_mD2->get_mD2(T);
                double lnQ2_1 = std::log(1. + taun / it->mfp0);
                double lnQ2_0 = std::log(1. + 
                    6.*it->p.boost_to(v3[0], v3[1], v3[2]).t()*T/mD2);
                double log_factor = std::sqrt(lnQ2_1/lnQ2_0);
                double LPM = std::min(it->mfp0 / taun * log_factor, 1.);
                // The final acceptance factor
                Acceptance = LPM * Running;

                if (Srandom::rejection(Srandom::gen) < Acceptance){   
                    if (it->pid==21) {
                        // mother parton id unchanged, semi-classical   
                        pIn.p = pIn.p - it->p;   
                        //f << pIn.x.x0() << " " << (pIn.p.t()+pIn.p.z())/2. << " " << (it->p.t()+it->p.z())/2. << std::endl;        
                    }
                    if (it->pid!=21) {
                        // mother parton id unchanged, 
                        // no semi-classical analog --> pair production

                        it->mass = pid2mass(it->pid);
                        pIn.mass = pid2mass(it->pid);
                        pIn.pid = -it->pid;
                        pIn.p = pIn.p - it->p;
                        it->p0 = it->p;                 
                        pIn.p.a[0] = std::sqrt(pIn.mass*pIn.mass+pIn.p.pabs2());
                        it->p.a[0] = std::sqrt(it->mass*it->mass+it->p.pabs2());
                        nonclassical = true;
                    }
                    pIn.p = put_on_shell(pIn.p, pIn.pid);
                    it->p = put_on_shell(it->p, it->pid);
                    pIn.p0 = pIn.p;
                    double Scale = std::sqrt(kt2n+mD2);
                    pIn.Q0 = Scale;
                    pIn.Q00 = Scale;
                    it->Q0 = Scale;
                    it->Q00 = Scale;
                    it->is_virtual = false;
                    pOut_list.push_back(*it);
                }
                // final remove it from the radlist
                it = pIn.radlist.erase(it);
            }
        }
        if (nonclassical) pIn.radlist.clear();
    }
    // Finally, remove those rad particle whose kinematic constrain is  
    // violated
    for (std::vector<particle>::iterator it = pIn.radlist.begin(); 
         it != pIn.radlist.end();){
        if(it->p.t()>=pIn.p.t()) it = pIn.radlist.erase(it);
        else it++;
    }


    // Add the mother parton to the end of the output list
    //also transform back its pid
    pOut_list.push_back(pIn);
    return pOut_list.size();
}


