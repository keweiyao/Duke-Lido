#include "Onium_Disso_dR.h"
#include "Onium_predefine.h"
#include "predefine.h"
#include <cmath>
#include <gsl/gsl_math.h>
#include "random.h"

//-----------Quarkonium gluon absorption-------------------
// QQbar + g --> Q + Qbar
// params[0] = v, params[1] = T, params[2] = M, params[3] = M
// v is quarkonium velocity w.r.t. hydrocell
// k1 = gamma*(1.+v)/T
// k2 = gamma*(1.-v)/T
// 1S State
double dRdq_1S_gluon(double q, void * params_){
    double * params = static_cast<double *>(params_);
    double v = params[0];
    double T = params[1];
    double M = params[2];
    double E1S = params[3];
    double a_B = params[4];
    double prel = sqrt(M*(q-E1S));
    double gamma = 1./std::sqrt(1. - v*v);
    // different technique for small v
    if (v > 1e-4){
        double k1 = gamma*(1.+v)/T; 
        double k2 = gamma*(1.-v)/T;
        
        return prefactor_disso_gluon * M * T /gamma/gamma/v
               * q*q*prel * (fac1(q*k1) - fac1(q*k2)) 
               * Matrix1S(prel,a_B);
    }
    else{ // small v formula
        return prefactor_disso_gluon_small_v * M /gamma
               * pow(q,3.)*prel * nB(q/T) 
               * Matrix1S(prel,a_B);
    }
}

// 2S State
double dRdq_2S_gluon(double q, void * params_){
    double * params = static_cast<double *>(params_);
    double v = params[0];
    double T = params[1];
    double M = params[2];
    double E2S = params[3];
    double a_B = params[4];
    double prel = sqrt(M*(q-E2S));
    double gamma = 1./std::sqrt(1. - v*v);
    // different technique for small v
    if (v > 1e-4){
        double k1 = gamma*(1.+v)/T; 
        double k2 = gamma*(1.-v)/T;
        
        return prefactor_disso_gluon * M * T /gamma/gamma/v
               * q*q*prel * (fac1(q*k1) - fac1(q*k2)) 
               * Matrix2S(prel, a_B);
    }
    else{ // small v formula
        return prefactor_disso_gluon_small_v * M / gamma
               * pow(q,3.)*prel * nB(q/T) 
               * Matrix2S(prel, a_B);
    }
}


// 1P State
double dRdq_1P_gluon(double q, void * params_){
    double * params = static_cast<double *>(params_);
    double v = params[0];
    double T = params[1];
    double M = params[2];
    double E1P = params[3];
    double a_B = params[4];
    double prel = sqrt(M*(q-E1P));
    double gamma = 1./std::sqrt(1. - v*v);
    // different technique for small v
    if (v > 1e-4){

        double k1 = gamma*(1.+v)/T; 
        double k2 = gamma*(1.-v)/T;
        
        return prefactor_disso_gluon * M * T /gamma/gamma/v
               * q*q*prel * (fac1(q*k1) - fac1(q*k2)) 
               * Matrix1P(prel, a_B);
    }
    else{ // small v formula
        return prefactor_disso_gluon_small_v * M /gamma
               * pow(q,3.)*prel * nB(q/T) 
               * Matrix1P(prel, a_B);
    }
}

//-----------Quarkonium inelastic breakup-------------------
// 1S State
// QQbar + q --> Q + Qbar + q
double dRdp1dp2_1S_decay_ineq(double x[5], std::size_t dim, void * params_){
    double p1 = x[0];
    double c1 = x[1];
    double Prel = x[2];
    double c2 = x[3];
    double phi = x[4];
    
    double * params = static_cast<double *>(params_);
    double v = params[0];
    double T = params[1];
    double M = params[2];
    double E1S = params[3];
    double a_B = params[4];
    double p2 = p1 - E1S - Prel*Prel/M;
    if (p2 <= 0.0){
        return 0.0;
    }
    else{
        double gamma = 1./std::sqrt(1.-v*v);
        double s1 = std::sqrt(1.-c1*c1);
        double s2 = std::sqrt(1.-c2*c2);
        double phase1 = p1*nF(gamma*p1*(1.+v*c1)/T);
        double phase2 = p2*nFminus1(gamma*p2*(1.+v*c2)/T);
        double part_angle = s1*s2*std::cos(phi)+c1*c2;
        double prop = p1*p2*(1.+part_angle)/( p1*p1+p2*p2-2.*p1*p2*part_angle );
        return phase1 * phase2 * Prel*Prel * Matrix1S(Prel, a_B) * prop;
        // omit a prefactor and a 1/gamma, add them in the rate calculation
    }
}
// 1S State
// QQbar + q --> Q + Qbar + g
double dRdq1dq2_1S_decay_ineg(double x[5], std::size_t dim, void * params_){
    // do a change of variable from p2 to Prel
    double q1 = x[0];
    double c1 = x[1];
    double Prel = x[2];
    double c2 = x[3];
    double phi = x[4];

    double * params = static_cast<double *>(params_);
    double v = params[0];
    double T = params[1];
    double M = params[2];
    double E1S = params[3];
    double a_B = params[4];
    double q2 = q1 - E1S - Prel*Prel/M;
    if (q2 <= 0.0){
        return 0.0;
    }
    else{
        double gamma = 1./std::sqrt(1.-v*v);
        double s1 = std::sqrt(1.-c1*c1);
        double s2 = std::sqrt(1.-c2*c2);
        double phase1 = q1*nB(gamma*q1*(1.+v*c1)/T);
        double phase2 = q2*nBplus1(gamma*q2*(1.+v*c2)/T);
        double part_angle = s1*s2*std::cos(phi)+c1*c2;
        double q1q2sum = q1+q2;
        double prop = q1q2sum*q1q2sum*(1.+part_angle)/( q1*q1+q2*q2-2.*q1*q2*part_angle );
        return phase1 * phase2 * Prel*Prel * Matrix1S(Prel, a_B) * prop;
        // omit a prefactor and a 1/gamma, add them in the rate calculation
    }
}

// 2S State
// QQbar + q --> Q + Qbar + q
double dRdp1dp2_2S_decay_ineq(double x[5], std::size_t dim, void * params_){
    double p1 = x[0];
    double c1 = x[1];
    double Prel = x[2];
    double c2 = x[3];
    double phi = x[4];

    double * params = static_cast<double *>(params_);
    double v = params[0];
    double T = params[1];
    double M = params[2];
    double E2S = params[3];
    double a_B = params[4];
    double p2 = p1 - E2S - Prel*Prel/M;
    if (p2 <= 0.0){
        return 0.0;
    }
    else{
        double gamma = 1./std::sqrt(1.-v*v);
        double s1 = std::sqrt(1.-c1*c1);
        double s2 = std::sqrt(1.-c2*c2);
        double phase1 = p1*nF(gamma*p1*(1.+v*c1)/T);
        double phase2 = p2*nFminus1(gamma*p2*(1.+v*c2)/T);
        double part_angle = s1*s2*std::cos(phi)+c1*c2;
        double prop = p1*p2*(1.+part_angle)/( p1*p1+p2*p2-2.*p1*p2*part_angle );
        return phase1 * phase2 * Prel*Prel * Matrix2S(Prel, a_B) * prop;
        // omit a prefactor and a 1/gamma, add them in the rate calculation
    }
}
// 2S State
// QQbar + g --> Q + Qbar + g
double dRdq1dq2_2S_decay_ineg(double x[5], std::size_t dim, void * params_){
    // do a change of variable from p2 to Prel
    double q1 = x[0];
    double c1 = x[1];
    double Prel = x[2];
    double c2 = x[3];
    double phi = x[4];

    double * params = static_cast<double *>(params_);
    double v = params[0];
    double T = params[1];
    double M = params[2];
    double E2S = params[3];
    double a_B = params[4];
    double q2 = q1 - E2S - Prel*Prel/M;
    if (q2 <= 0.0){
        return 0.0;
    }
    else{
        double gamma = 1./std::sqrt(1.-v*v);
        double s1 = std::sqrt(1.-c1*c1);
        double s2 = std::sqrt(1.-c2*c2);
        double phase1 = q1*nB(gamma*q1*(1.+v*c1)/T);
        double phase2 = q2*nBplus1(gamma*q2*(1.+v*c2)/T);
        double part_angle = s1*s2*std::cos(phi)+c1*c2;
        double q1q2sum = q1+q2;
        double prop = q1q2sum*q1q2sum*(1.+part_angle)/( q1*q1+q2*q2-2.*q1*q2*part_angle );
        return phase1 * phase2 * Prel*Prel * Matrix2S(Prel, a_B) * prop;
        // omit a prefactor and a 1/gamma, add them in the rate calculation
    }
}


// 1P State
// QQbar + q --> Q + Qbar + q
double dRdp1dp2_1P_decay_ineq(double x[5], std::size_t dim, void * params_){
    double p1 = x[0];
    double c1 = x[1];
    double Prel = x[2];
    double c2 = x[3];
    double phi = x[4];

    double * params = static_cast<double *>(params_);
    double v = params[0];
    double T = params[1];
    double M = params[2];
    double E1P = params[3];
    double a_B = params[4];
    double p2 = p1 - E1P - Prel*Prel/M;
    if (p2 <= 0.0){
        return 0.0;
    }
    else{
        double gamma = 1./std::sqrt(1.-v*v);
        double s1 = std::sqrt(1.-c1*c1);
        double s2 = std::sqrt(1.-c2*c2);
        double phase1 = p1*nF(gamma*p1*(1.+v*c1)/T);
        double phase2 = p2*nFminus1(gamma*p2*(1.+v*c2)/T);
        double part_angle = s1*s2*std::cos(phi)+c1*c2;
        double prop = p1*p2*(1.+part_angle)/( p1*p1+p2*p2-2.*p1*p2*part_angle );
        return phase1 * phase2 * Prel*Prel * Matrix1P(Prel, a_B) * prop;
        // omit a prefactor and a 1/gamma, add them in the rate calculation
    }
}
// 1P State
// QQbar + g --> Q + Qbar + g
double dRdq1dq2_1P_decay_ineg(double x[5], std::size_t dim, void * params_){
    // do a change of variable from p2 to Prel
    double q1 = x[0];
    double c1 = x[1];
    double Prel = x[2];
    double c2 = x[3];
    double phi = x[4];

    double * params = static_cast<double *>(params_);
    double v = params[0];
    double T = params[1];
    double M = params[2];
    double E1P = params[3];
    double a_B = params[4];
    double q2 = q1 - E1P - Prel*Prel/M;
    if (q2 <= 0.0){
        return 0.0;
    }
    else{
        double gamma = 1./std::sqrt(1.-v*v);
        double s1 = std::sqrt(1.-c1*c1);
        double s2 = std::sqrt(1.-c2*c2);
        double phase1 = q1*nB(gamma*q1*(1.+v*c1)/T);
        double phase2 = q2*nBplus1(gamma*q2*(1.+v*c2)/T);
        double part_angle = s1*s2*std::cos(phi)+c1*c2;
        double q1q2sum = q1+q2;
        double prop = q1q2sum*q1q2sum*(1.+part_angle)/( q1*q1+q2*q2-2.*q1*q2*part_angle );
        return phase1 * phase2 * Prel*Prel * Matrix1P(Prel, a_B) * prop;
        // omit a prefactor and a 1/gamma, add them in the rate calculation
    }
}

// -------------------------- Sampling --------------------------
// disso_real_gluon
double disso_gluon_costheta(double q, double v, double T, double random){
    v = std::max(v, 1e-4);
    double coeff = q/T/std::sqrt(1. - v*v);     // parameter B = coeff
    double low = fac1(coeff*(1.-v));
    double norm = fac1(coeff*(1.+v))-low;
    double y_fac1 = random*norm + low;
    return -(1. + std::log(1. - std::exp(y_fac1))/coeff )/v;
}

// inelastic quark scattering disso
 // function used in importance sampling of p1
 double f_p1_disso_important(double p1, void * params_){
     double * params = static_cast<double *>(params_);
     double v = params[0];
     double T = params[1];
     double Enl = params[2];
     double gamma = 1./std::sqrt(1.-v*v);
     double k1 = gamma*(1-v)/T;        // k1 = gamma*(1-v)/T
     double k2 = gamma*(1+v)/T;        // k2 = gamma*(1+v)/T
     double p2 = p1-Enl;     // p2 = p1 - Enl
     if (p2 <= 0){
         return 0.0;
     }
     else{
         double p1p2 = p1/p2;
         return ( fac2(k1*p1) - fac2(k2*p1) ) * p1 /(p1p2 + 1./p1p2 -2.0);
     }
 }

double Sample_disso_ineq_p1_important(double p1low, double p1up, double result_max, void * params_){ // result_max is an input
     double * params = static_cast<double *>(params_);
     double result_try, p1_try;
     do{
         p1_try = Srandom::rejection(Srandom::gen)*(p1up-p1low) + p1low;
         result_try = f_p1_disso_important(p1_try, params);
     } while(Srandom::rejection(Srandom::gen)*result_max > result_try);
     return p1_try;
}
 
double Sample_disso_ineq_cos1(double p1, void * params_){
    double * params = static_cast<double *>(params_);
    double y_try = Srandom::y_cdf(Srandom::gen);
    double v = std::max(params[0], 1e-4);
    double B = params[1]*p1;    // B = gamma * p1/T
    double C = y_try * fac2(B*(1.+v)) + (1.-y_try) * fac2(B*(1.-v));
    return -(1. + std::log(std::exp(C) - 1.)/B )/v;
 }
 
 std::vector<double> Sample_disso_ineq(double v, double T, double mass, double Enl, double a_B, double prel_up, double maximum_f_p1, double max_p2Matrix, double(*f)(double _prel, double _aB)){    // maximum is the input for Sample_disso_ineq_p1_important
     v = std::max(v, 1e-4);
     double gamma = 1./std::sqrt(1.-v*v);
     double p1_try, c1_try, s1_try, p2_try, c2_try, s2_try, phi_try, c_phi, s_phi, p_rel, p1p2, part_angle, result_try, p1p2_try;
     double p1low = Enl;
     double p1up = 15.*T/std::sqrt(1.-v);
 
     double * params_p1 = new double[3];
     params_p1[0] = v; //k1 = gamma*(1-v)/T
     params_p1[1] = T; //k2 = gamma*(1+v)/T
     params_p1[2] = Enl;
 
     double * params_c1 = new double[2];
     params_c1[0] = v;
     params_c1[1] = gamma/T;
 
     do{
         do{
             p_rel = Srandom::sample_inel(Srandom::gen)*prel_up;
         } while(Srandom::rejection(Srandom::gen)*max_p2Matrix > f(p_rel, a_B));
 
         p1_try = Sample_disso_ineq_p1_important(p1low, p1up, maximum_f_p1, params_p1);  // give the maximum as result_max to Sample_disso_ineq_p1_important
         p2_try = p1_try - Enl - p_rel*p_rel/mass;
         if (p2_try <= 0.0){
             result_try = 0.0;
         }
         else{
             c1_try = Sample_disso_ineq_cos1(p1_try, params_c1);
             c2_try = Srandom::dist_costheta(Srandom::gen);
             s1_try = std::sqrt(1.-c1_try*c1_try);
             s2_try = std::sqrt(1.-c2_try*c2_try);
             phi_try = Srandom::dist_phi(Srandom::gen);
             p1p2_try = p1_try/p2_try;   // p1_try/p2_try
             p1p2 = (p1_try-Enl)/p1_try;
             c_phi = std::cos(phi_try);
             s_phi = std::sin(phi_try);
             part_angle = s1_try*s2_try*c_phi + c1_try*c2_try;
             result_try = (1.+part_angle)/(p1p2_try + 1./p1p2_try - 2.*part_angle)/2.0*(p1p2 + 1./p1p2 - 2.0) * nFminus1(gamma*(1.+v*c2_try)*p2_try/T)/p1p2_try;
         }
     } while(Srandom::rejection(Srandom::gen) >= result_try);
     std::vector<double> p1_final(3);
     std::vector<double> p2_final(3);
     std::vector<double> p_rel_final(3);
     std::vector<double> pQpQbar_final(6);
     double cos_rel, phi_rel;
     cos_rel = Srandom::dist_costheta(Srandom::gen);
     phi_rel = Srandom::sample_inel(Srandom::gen)*TwoPi;
     p1_final = polar_to_cartisian2(p1_try, c1_try, s1_try, 1.0, 0.0);
     p2_final = polar_to_cartisian2(p2_try, c2_try, s2_try, c_phi, s_phi);
     p_rel_final = polar_to_cartisian1(p_rel, cos_rel, phi_rel);
     pQpQbar_final = add_virtual_gluon(p1_final, p2_final, p_rel_final);
     return pQpQbar_final;
 }

// inelastic gluon scattering disso
    // function used in importance sampling of q1
double f_q1_disso_important(double q1, void * params_){
    double * params = static_cast<double *>(params_);
    double v = params[0];
    double T = params[1];
    double Enl = params[2];
    double gamma = 1./std::sqrt(1.-v*v);
    double k1 = gamma*(1-v)/T;        // k1 = gamma*(1-v)/T
    double k2 = gamma*(1+v)/T;        // k2 = gamma*(1+v)/T
    double q2 = q1-Enl;     // p2 = p1 - Enl
    if (q2 <= 0){
        return 0.0;
    }
    else{
        double q1q2 = q1/q2;
        return ( fac1(k2*q1) - fac1(k1*q1) ) * q1 * ( 1.0 + 4./(q1q2 + 1./q1q2 -2.0) );
    }
}

double Sample_disso_ineg_q1_important(double q1low, double q1up, double result_max, void * params_){ // result_max is an input
    double * params = static_cast<double *>(params_);
    double result_try, q1_try;
    do{
        q1_try = Srandom::rejection(Srandom::gen)*(q1up-q1low) + q1low;
        result_try = f_q1_disso_important(q1_try, params);
    } while(Srandom::rejection(Srandom::gen)*result_max > result_try);
    return q1_try;
}

double Sample_disso_ineg_cos1(double q1, void * params_){
    double * params = static_cast<double *>(params_);
    double y_try = Srandom::y_cdf(Srandom::gen);
    double v = std::max(params[0], 1e-4);
    double B = params[1]*q1;    // B = gamma * q1/T
    double C = y_try * fac1(B*(1.+v)) + (1.-y_try) * fac1(B*(1.-v));
    return -(1. + std::log( 1.-std::exp(C) )/B )/v;
}

std::vector<double> Sample_disso_ineg(double v, double T, double mass, double Enl, double a_B, double prel_up, double maximum_f_q1, double max_p2Matrix, double(*f)(double _prel, double _aB)){    // maximum is the input for Sample_disso_ineq_q1_important
    v = std::max(v, 1e-4);
    double gamma = 1./std::sqrt(1.-v*v);
    double q1_try, c1_try, s1_try, q2_try, c2_try, s2_try, phi_try, c_phi, s_phi, p_rel, q1q2, q1q2_try, part_angle, result_try;
    double q1low = Enl;
    double q1up = 15.*T/std::sqrt(1.-v);
    
    double * params_q1 = new double[3];
    params_q1[0] = v; //k1 = gamma*(1-v)/T
    params_q1[1] = T; //k2 = gamma*(1+v)/T
    params_q1[2] = Enl;
    
    double * params_c1 = new double[2];
    params_c1[0] = v;
    params_c1[1] = gamma/T;
    
    do{
        do{
            p_rel = Srandom::sample_inel(Srandom::gen)*prel_up;
        } while(Srandom::rejection(Srandom::gen)*max_p2Matrix > f(p_rel, a_B));
        
        q1_try = Sample_disso_ineg_q1_important(q1low, q1up, maximum_f_q1, params_q1);  // give the maximum as result_max to Sample_disso_ineg_q1_important
        q2_try = q1_try - Enl - p_rel*p_rel/mass;
        if (q2_try <= 0.0){
            result_try = 0.0;
        }
        else{
            c1_try = Sample_disso_ineg_cos1(q1_try, params_c1);
            c2_try = Srandom::dist_costheta(Srandom::gen);
            s1_try = std::sqrt(1.-c1_try*c1_try);
            s2_try = std::sqrt(1.-c2_try*c2_try);
            phi_try = Srandom::dist_phi(Srandom::gen);
            q1q2_try = q1_try/q2_try;   // q1_try/q2_try
            q1q2 = (q1_try-Enl)/q1_try;
            c_phi = std::cos(phi_try);
            s_phi = std::sin(phi_try);
            part_angle = s1_try*s2_try*c_phi + c1_try*c2_try;
            result_try = nBplus1(gamma*(1.+v*c2_try)*q2_try/T)/q1q2_try * (q1q2_try  + 1./q1q2_try + 2.0)/(q1q2_try + 1./q1q2_try - 2.0*part_angle) * (1.+part_angle)/2.0 / (1.0 + 4./(q1q2 + 1./q1q2 - 2.0) );
        }
    } while(Srandom::rejection(Srandom::gen) >= result_try);
    std::vector<double> q1_final(3);
    std::vector<double> q2_final(3);
    std::vector<double> p_rel_final(3);
    std::vector<double> pQpQbar_final(6);
    double cos_rel, phi_rel;
    cos_rel = Srandom::dist_costheta(Srandom::gen);
    phi_rel = Srandom::sample_inel(Srandom::gen)*TwoPi;
    q1_final = polar_to_cartisian2(q1_try, c1_try, s1_try, 1.0, 0.0);
    q2_final = polar_to_cartisian2(q2_try, c2_try, s2_try, c_phi, s_phi);
    p_rel_final = polar_to_cartisian1(p_rel, cos_rel, phi_rel);
    pQpQbar_final = add_virtual_gluon(q1_final, q2_final, p_rel_final);
    return pQpQbar_final;
}
