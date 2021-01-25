#include "collision_manager.h"
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/variant/get.hpp>
#include <boost/any.hpp>
#include <boost/foreach.hpp>
#include "matrix_elements.h"
#include "simpleLogger.h"
#include "Langevin.h"
#include "random.h"

collision_manager::collision_manager(std::string mode,
                    std::string setting_path,
                    std::string table_path, 
                    std::vector<double> parameters):
AllProcesses(){
    double mu = parameters[0], afix = parameters[1], 
           cut = parameters[2], theta = parameters[3];
    initialize_mD_and_scale(0, mu, afix, theta, cut);
    echo();
   // For gluon, 17 channels
    AllProcesses[21] = std::vector<Process>();
    AllProcesses[21].push_back(Rate22("Boltzmann/gq2gq", 
                               setting_path, dX_gq2gq));    
    AllProcesses[21].push_back(Rate23("Boltzmann/gq2gqg", 
                               setting_path, dX_gq2gqg));  
    AllProcesses[21].push_back(Rate23("Boltzmann/gq2qqqbar",
                               setting_path, dX_gq2qqqbar)); 
    AllProcesses[21].push_back(Rate23("Boltzmann/gq2cqcbar", 
                               setting_path, dX_gq2qqqbar)); 
    AllProcesses[21].push_back(Rate23("Boltzmann/gq2bqbbar", 
                               setting_path, dX_gq2qqqbar)); 
    AllProcesses[21].push_back(Rate22("Boltzmann/gg2gg", 
                               setting_path, dX_gg2gg));      
    AllProcesses[21].push_back(Rate22("Boltzmann/gg2qqbar", 
                               setting_path, dX_gg2qqbar));   
    AllProcesses[21].push_back(Rate22("Boltzmann/gg2ccbar", 
                               setting_path, dX_gg2qqbar));   
    AllProcesses[21].push_back(Rate22("Boltzmann/gg2bbbar", 
                               setting_path, dX_gg2qqbar));   
    AllProcesses[21].push_back(Rate23("Boltzmann/gg2ggg", 
                               setting_path, dX_gg2ggg));      
    AllProcesses[21].push_back(Rate23("Boltzmann/gg2qgqbar", 
                               setting_path, dX_gg2qgqbar));   
    AllProcesses[21].push_back(Rate23("Boltzmann/gg2cgcbar", 
                               setting_path, dX_gg2qgqbar)); 
    AllProcesses[21].push_back(Rate23("Boltzmann/gg2bgbbar", 
                               setting_path, dX_gg2qgqbar)); 
    AllProcesses[21].push_back(Rate23("Boltzmann/gg2gqqbar", 
                               setting_path, dX_gg2gqqbar));   
    AllProcesses[21].push_back(Rate23("Boltzmann/gg2gccbar", 
                               setting_path, dX_gg2gqqbar)); 
    AllProcesses[21].push_back(Rate23("Boltzmann/gg2gbbbar", 
                               setting_path, dX_gg2gqqbar));     
    AllProcesses[21].push_back(Rate12("Boltzmann/g2gg", 
                               setting_path, LGV_g2gg));       
    AllProcesses[21].push_back(Rate12("Boltzmann/g2qqbar", 
                               setting_path, LGV_g2qqbar));  
    AllProcesses[21].push_back(Rate12("Boltzmann/g2ccbar", 
                               setting_path, LGV_g2qqbar));   
    AllProcesses[21].push_back(Rate12("Boltzmann/g2bbbar", 
                               setting_path, LGV_g2qqbar)); 
    // For light flavor, 11 channels
    AllProcesses[123] = std::vector<Process>();
    AllProcesses[123].push_back(Rate22("Boltzmann/qq2qq", 
                                setting_path, dX_qq2qq)); 
    AllProcesses[123].push_back(Rate23("Boltzmann/qq2qqg", 
                                setting_path, dX_qq2qqg));  
    AllProcesses[123].push_back(Rate22("Boltzmann/qg2qg", 
                                setting_path, dX_qg2qg));
    AllProcesses[123].push_back(Rate23("Boltzmann/qg2qgg", 
                                setting_path, dX_qg2qgg));  
    AllProcesses[123].push_back(Rate23("Boltzmann/qg2qqqbar", 
                                setting_path, dX_qg2qqqbar)); 
    AllProcesses[123].push_back(Rate23("Boltzmann/qg2qccbar", 
                                setting_path, dX_qg2qqqbar)); 
    AllProcesses[123].push_back(Rate23("Boltzmann/qg2qbbbar", 
                                setting_path, dX_qg2qqqbar)); 
    AllProcesses[123].push_back(Rate22("Boltzmann/qqbar2qqbar", 
                                setting_path, dX_qqbar2qqbar_diff)); 
    AllProcesses[123].push_back(Rate22("Boltzmann/qqbar2ccbar", 
                                setting_path, dX_qqbar2qqbar_diff)); 
    AllProcesses[123].push_back(Rate22("Boltzmann/qqbar2bbbar", 
                                setting_path, dX_qqbar2qqbar_diff)); 
    AllProcesses[123].push_back(Rate12("Boltzmann/q2qg", 
                                setting_path, LGV_q2qg));     
    // for open charm flavor 5 channels
    AllProcesses[4] = std::vector<Process>();
    AllProcesses[4].push_back(Rate22("Boltzmann/cq2cq", 
                              setting_path, dX_qq2qq)); 
    AllProcesses[4].push_back(Rate22("Boltzmann/cg2cg", 
                              setting_path, dX_qg2qg)); 
    AllProcesses[4].push_back(Rate23("Boltzmann/cq2cqg", 
                              setting_path, dX_qq2qqg)); 
    AllProcesses[4].push_back(Rate23("Boltzmann/cg2cgg", 
                              setting_path, dX_qg2qgg)); 
    AllProcesses[4].push_back(Rate12("Boltzmann/c2cg", 
                              setting_path, LGV_q2qg));  
    // for open bottom flavor 5 channels
    AllProcesses[5] = std::vector<Process>();
    AllProcesses[5].push_back(Rate22("Boltzmann/bq2bq", 
                              setting_path, dX_qq2qq)); 
    AllProcesses[5].push_back(Rate22("Boltzmann/bg2bg", 
                              setting_path, dX_qg2qg));
    AllProcesses[5].push_back(Rate23("Boltzmann/bq2bqg", 
                              setting_path, dX_qq2qqg)); 
    AllProcesses[5].push_back(Rate23("Boltzmann/bg2bgg", 
                              setting_path, dX_qg2qgg)); 
    AllProcesses[5].push_back(Rate12("Boltzmann/b2bg", 
                              setting_path, LGV_q2qg));

    // Initialzie all processes
    std::vector<int> IDS{21,123,4,5};
    for (auto & i : IDS){
        BOOST_FOREACH(Process &r, AllProcesses[i])
            init_process(r, mode, table_path);
    }
}

double collision_manager::get_effective_mfp_soft(
                              double Ecell, double x,
                              int idA, int idB, int idC,
                              double Temp
                              ){
    double qA = qhat(idA, Ecell, 0., Temp);
    double qB = qhat(idB, x*Ecell, 0., Temp);
    double qC = qhat(idC, (1-x)*Ecell, 0., Temp);
    double mD2 = t_channel_mD2->get_mD2(Temp);
    double qhat3 = x*qC+(1.-x)*qB-x*(1.-x)*qA;
    return LPM_prefactor * mD2 /qhat3;
}

double collision_manager::get_effective_mfp_hard(
                              double Ecell, double x,
                              int idA, int idB, int idC,
                              double Temp){
    //double qA=0., qB=0., qC=0.;
    double qA = qhat(idA, Ecell, 0., Temp);
    double qB = qhat(idB, x*Ecell, 0., Temp);
    double qC = qhat(idC, (1-x)*Ecell, 0., Temp);
    double qhat3soft = x*qC+(1.-x)*qB-x*(1.-x)*qA;
    double mD2 = t_channel_mD2->get_mD2(Temp);
    double lnEcell = std::log(Ecell);

    int iA = std::abs(idA),iB = std::abs(idB),iC = std::abs(idC);
    if (iA<=3) iA = 123;
    if (iB<=3) iB = 123;
    if (iC<=3) iC = 123;
    double CRA = (idA==21)? CA:CF;
    double CRB = (idB==21)? CA:CF;
    double CRC = (idC==21)? CA:CF;

    BOOST_FOREACH (Process &r, AllProcesses[iA]){
        switch (r.which()){
        case 0:
            if (boost::get<Rate22>(r).IsActive()){
                int prcid = boost::get<Rate22>(r).process_id();
                if (is_classical(prcid)) {
                    auto A = boost::get<Rate22>(r).GetSecondM(
                               {std::log(Ecell), Temp});
                    qA += A.T[1][1] + A.T[2][2];
                }
            }
            break;
        default:
            break;
        }
    }  
    BOOST_FOREACH (Process &r, AllProcesses[iB]){
        switch (r.which()){
        case 0:
            if (boost::get<Rate22>(r).IsActive()){
                int prcid = boost::get<Rate22>(r).process_id();
                if (is_classical(prcid)) {
                    auto A = boost::get<Rate22>(r).GetSecondM(
                               {std::log(x*Ecell), Temp});
                    qB += A.T[1][1] + A.T[2][2];
                }
            }
            break;
        default:
            break;
        }
    }  
    BOOST_FOREACH (Process &r, AllProcesses[iC]){
        switch (r.which()){
        case 0:
            if (boost::get<Rate22>(r).IsActive()){
                int prcid = boost::get<Rate22>(r).process_id();
                if (is_classical(prcid)) {
                    auto A = boost::get<Rate22>(r).GetSecondM(
                               {std::log((1-x)*Ecell), Temp});
                    qC += A.T[1][1] + A.T[2][2];
                }
            }
            break;
        default:
            break;
        }
    }        
    double qhat3hard = x*qC+(1.-x)*qB-x*(1.-x)*qA;
    return LPM_prefactor * mD2/ (qhat3hard+qhat3soft) 
           * (x*CRC+(1-x)*CRB-x*(1-x)*CRA) /CA
                               ;
}

void collision_manager::init_process(Process &r, std::string mode, std::string table_path){
    switch (r.which()){
    case 0:
        if (boost::get<Rate22>(r).IsActive())
            //LOG_INFO << boost::get<Rate22>(r).name
            if (mode == "new"){
                boost::get<Rate22>(r).initX(table_path);
                boost::get<Rate22>(r).init(table_path);
            }
            else{
                boost::get<Rate22>(r).loadX(table_path);
                boost::get<Rate22>(r).load(table_path);
            }
        break;
    case 1:
        if (boost::get<Rate23>(r).IsActive())
            if (mode == "new"){
                boost::get<Rate23>(r).initX(table_path);
                boost::get<Rate23>(r).init(table_path);
            }
            else{
                boost::get<Rate23>(r).loadX(table_path);
                boost::get<Rate23>(r).load(table_path);
            }
        break;
    case 2:
        if (boost::get<Rate12>(r).IsActive()){
            if (mode == "new") boost::get<Rate12>(r).init(table_path);
            else boost::get<Rate12>(r).load(table_path);
        }
        break;
    default:
        LOG_FATAL << "Initializing processes: process type" 
                  << r.which() << " not exists";
        exit(-1);
        break;
    }
}

int collision_manager::sample(particle & pIn, 
                              std::vector<particle> & pOut_list,
                              double Temp, 
                              std::vector<double> v3cell,
                              double dt_lab){
    fourvec p0 = pIn.p;
    pOut_list.clear();
    // check pid
    int HardID;
    int abspid = std::abs(pIn.pid);
    if (abspid == 21) HardID = 21;
    else if (abspid==1 || abspid==2 || abspid==3) HardID = 123;
    else if (abspid==4) HardID = 4;
    else if (abspid==5) HardID = 5;
    else {
        LOG_FATAL << "1. This particle, pid="
                  << pIn.pid
                  << " should not be handled by the LIDO collision term!";
        exit(-1);
    }    

    // Get the kinematics in the cell frame:
    auto p_cell = pIn.p.boost_to(v3cell[0], v3cell[1], v3cell[2]);
    double E_cell = std::max(p_cell.t(), .5);
    double LnEcell = std::log(E_cell);
    double dt_cell = dt_lab / pIn.p.t() * p_cell.t();
    // For diffusion induced radiation, qhat_Soft is an input
    double softqhatg = qhat(21, E_cell, 0., Temp);

    // compute total incoherent collision rate in the cell frame
    std::vector<double> CumProb(AllProcesses[HardID].size());
    double TotProb = 0., dR;
    bool IsReal = (!pIn.is_virtual); // only real particles has all processes
                                     // virtual ones only do classical ones
                                     // : t-channel, elastic, ID preserving
    int count = 0;
    BOOST_FOREACH (Process &r, AllProcesses[HardID]){
        switch (r.which()){
        case 0:
            if (boost::get<Rate22>(r).IsActive())
                if (IsReal)
                    dR = boost::get<Rate22>(r).GetZeroM({LnEcell, Temp}).s;
                else{
                    int prcid = boost::get<Rate22>(r).process_id();
                    if (is_classical(prcid))
                        dR = boost::get<Rate22>(r).GetZeroM({LnEcell, Temp}).s;
                    else dR = 0.;
                }
            else dR = 0.0;
            CumProb[count] = TotProb + dR * dt_cell;
            break;
        case 1:
            if (boost::get<Rate23>(r).IsActive() && IsReal)
                dR = boost::get<Rate23>(r).GetZeroM({LnEcell, Temp}).s;
            else dR = 0.0;
            CumProb[count] = TotProb + dR * dt_cell;
            break;
        case 2:
            if (boost::get<Rate12>(r).IsActive() && IsReal)
                dR = softqhatg*boost::get<Rate12>(r).GetZeroM({LnEcell, Temp}).s;
            else dR = 0.0;
            CumProb[count] = TotProb + dR * dt_cell;
            break;
        default:
            LOG_FATAL << "Computing rate: process type" 
                      << r.which() << " not exists";
            exit(-1);
            break;
        }
        TotProb += dR*dt_cell;
        count++;
    }
    // P_total is the total rate of reation within dt
    if (TotProb > 0.15 && !type2_warned){
        LOG_WARNING << "P(dt) = " << TotProb << " may be too large";
        type2_warned = true;
    }
    // sample collision channel
    int choice = -1;
    if (Srandom::init_dis(Srandom::gen) > TotProb) {
        pOut_list.push_back(pIn);   
        return -1; 
    }
    else{
        std::vector<double> NormProb(AllProcesses[HardID].size());
        // Normalize P_channels using P_total to sample
        for (int i = 0; i < CumProb.size(); ++i) NormProb[i] = CumProb[i]/TotProb;
        double p = Srandom::init_dis(Srandom::gen);
        for (int i = 0; i < NormProb.size(); ++i){
            if (NormProb[i] > p){
                choice = i;
                break;
            }
        }
    }

    // Sample incoherent final state in the cell frame
    if (choice >= AllProcesses[HardID].size()){
        LOG_INFO << pIn.pid << " p0 = " << p0 << ", p=" << pIn.p;
        LOG_INFO << Temp  << " " 
                << v3cell[0] << " " << v3cell[1] << " " << v3cell[2];
        LOG_FATAL << "Choosen channel: " << choice << " does not exist";
        exit(-1);
    }
    // Final state holder FS, and pids
    std::vector<fourvec> FS;
    std::vector<int> FS_pids;
    // Sampe differential final state
    int prcid = -3;
    bool status;
    switch (AllProcesses[HardID][choice].which()){
    case 0:
        status = boost::get<Rate22>(AllProcesses[HardID][choice]).sample(
                          {LnEcell, Temp}, pIn.pid, FS, FS_pids);
        prcid = boost::get<Rate22>(AllProcesses[HardID][choice]).process_id();
        break;
    case 1:
        status = boost::get<Rate23>(AllProcesses[HardID][choice]).sample(
                          {LnEcell, Temp}, pIn.pid, FS, FS_pids);
        prcid = boost::get<Rate23>(AllProcesses[HardID][choice]).process_id();
        break;
    case 2:
        status = boost::get<Rate12>(AllProcesses[HardID][choice]).sample(
                          {LnEcell, Temp}, pIn.pid, FS, FS_pids);
        prcid = boost::get<Rate12>(AllProcesses[HardID][choice]).process_id();
        break;
    default:
        LOG_FATAL << "Sampling rate: process type" 
                  << choice << " not exists";
        exit(-1);
        break;
    }
    // status==false means sampling has failed,
    if (!status){
        LOG_INFO << "Collision manager: sampling failed";
        pOut_list.clear();   
        return -2;    
    }
    // if the particle is virtual, rescale its energy
    if (pIn.is_virtual) {
        FS[0] = FS[0]*(p_cell.t()/FS[0].t());
        FS[0] = put_on_shell(FS[0], FS_pids[0]);
    }
    // rotate the final state back and boost it back to the lab frame
    for (auto &pmu : FS){
        pmu = pmu.rotate_back(p_cell);
        pmu = pmu.boost_back(v3cell[0], v3cell[1], v3cell[2]);
    }
    // organize final states
    pOut_list.clear();
    // if it is a 2->2 process
    if(is2to2(prcid)){
        for (int i=0; i<FS.size(); i++){
            int col=0, anticol=0;
            auto newp = make_parton(FS_pids[i], col, anticol, pIn.Q0,
                                FS[i], pIn, Temp, v3cell);
            pOut_list.push_back(newp);  
        }
        return prcid;
    }
    // if it is a 1->2 process
    double boost_factor = pIn.p.t() / E_cell;
    if(is1to2(prcid)){
        double p0abs = pIn.p.pabs();
        fourvec nbar{1., -pIn.p.x()/p0abs,  
               -pIn.p.y()/p0abs, -pIn.p.z()/p0abs};
        double x = dot(FS[1], nbar)/dot(pIn.p, nbar);
        double mfp = boost_factor * get_effective_mfp_hard(
              E_cell, x, pIn.pid, FS_pids[1], FS_pids[0], Temp);
        for (int i=0; i<FS.size(); i++){
            int col=0, anticol=0;
            auto newp = make_parton(FS_pids[i], col, anticol, pIn.Q0,
                                FS[i], pIn, Temp, v3cell);
	    newp.is_virtual = true;
            newp.mfp0 = mfp;
            pOut_list.push_back(newp);  
        }
        return prcid;
    }
    // if it is a 2->3 process
    if(is2to3(prcid)){
        double p0abs = pIn.p.pabs();
        fourvec nbar{1., -pIn.p.x()/p0abs,  
               -pIn.p.y()/p0abs, -pIn.p.z()/p0abs};
        double x = dot(FS[2], nbar)/dot(pIn.p, nbar);
        double mfp = boost_factor * get_effective_mfp_hard(
              E_cell, x, pIn.pid, FS_pids[2], FS_pids[0], Temp);
        for (int i=0; i<FS.size(); i++){
            int col=0, anticol=0;
            auto newp = make_parton(FS_pids[i], col, anticol, pIn.Q0,
                                FS[i], pIn, Temp, v3cell);
	    newp.is_virtual = true;
            newp.mfp0 = mfp;
            pOut_list.push_back(newp);  
        }
        return prcid;
    }
}

