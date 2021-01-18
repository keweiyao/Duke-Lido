#include "Langevin.h"
#include "random.h"
#include "predefine.h"
double const tiny = 1e-10;

double qhat_small_angle_LOpQCD(int pid, double E, double M, double T){
    double CR = (pid==21) ? CA : CF;
    double mD2 = t_channel_mD2->get_mD2(T);
    double Q2cut = std::max(std::min(cut*mD2, 6*E*T),mD2);
    return alpha_s(Q2cut, T) * CR * T * mD2 * std::log(1+Q2cut/mD2);
}

double qhat_L_small_angle_LOpQCD(int pid, double E, double M, double T){
    double CR = (pid==21) ? CA : CF;
    double mD2 = t_channel_mD2->get_mD2(T);
    double Q2cut = std::max(std::min(cut*mD2, 6*E*T),mD2);
    return alpha_s(Q2cut, T) * CR * T * .5*mD2 * std::log(1+Q2cut/mD2);
}


double qhat(int pid, double E, double M, double T){
    return  qhat_small_angle_LOpQCD(pid, E, M, T);
}


double qhat_L(int pid, double E, double M, double T){
    return  qhat_L_small_angle_LOpQCD(pid, E, M, T);
}

double dqhat_L_dp2(int pid, double E, double M, double T){
    double p2 = E*E - M*M + tiny;
    double dp2 = p2*.05;
    double Eprime = std::sqrt(E*E + dp2);
    return (qhat_L(pid, Eprime, M, T) - qhat_L(pid, E, M, T) ) /dp2;
}

void Ito_update_lab(int pid, double dt_lab, double M, double T, 
    std::vector<double> v, const fourvec & pIn, fourvec & pOut){
    // Boost pIn to medium frame
    auto pIn_cell = pIn.boost_to(v[0], v[1], v[2]);
    // Boost dt to medium frame
    double dt = dt_lab*pIn_cell.t()/pIn.t();
    // imaging rotating to a frame where pIn lies on z-axis
    double E0 = pIn_cell.t();
    double p0 = std::sqrt(E0*E0 - M*M + 1e-6);

    double kt = qhat(pid, E0, M, T)/2.;
    double kl = qhat_L(pid, E0, M, T);
    double Ed = std::max(E0, 3.*T);
    double drag = kl/(2.*Ed*T);
    double damp = std::max(1.-drag*dt, 0.);           
    double Ct = std::sqrt(kt*dt);
    double Cl = std::sqrt(kl*dt);
   
    pOut.a[1] = Ct * Srandom::white_noise(Srandom::gen);
    pOut.a[2] = Ct * Srandom::white_noise(Srandom::gen);
    pOut.a[3] = p0 * damp + Cl * Srandom::white_noise(Srandom::gen);
    pOut.a[0] = std::sqrt(M*M + pOut.pabs2() );

    // rotate back to the original frame
    pOut = pOut.rotate_back(pIn_cell);
    // boost back to lab frame
    pOut = pOut.boost_back(v[0], v[1], v[2]);
}

void Ito_update_rest(int pid, double dt, double M, double T, 
                     const fourvec & pIn, fourvec & pOut){
    // imaging rotating to a frame where pIn lies on z-axis
    double E0 = pIn.t();
    double p0 = std::sqrt(E0*E0 - M*M + 1e-9);
    double Eregulate = std::max(E0, 3*T);
    double kt = qhat(pid, Eregulate, M, T)/2.;
    double kl = qhat_L(pid, Eregulate, M, T);
    double drag = kl/(2.*Eregulate*T);
    double damp = std::max(1.-drag*dt, 0.);           
    double Ct = std::sqrt(kt*dt);
    double Cl = std::sqrt(kl*dt);
    pOut.a[1] = Ct * Srandom::white_noise(Srandom::gen);
    pOut.a[2] = Ct * Srandom::white_noise(Srandom::gen);
    pOut.a[3] = p0 * damp + Cl * Srandom::white_noise(Srandom::gen);
    pOut.a[0] = std::sqrt(M*M + pOut.pabs2() );
    // rotate back to the original frame
    pOut = pOut.rotate_back(pIn);
}

