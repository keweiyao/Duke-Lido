#include "Onium_Disso_dR.h"
#include "Onium_predefine.h"
#include "predefine.h"
#include <cmath>
#include <gsl/gsl_math.h>

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

