#include "Onium_Reco_dR.h"
#include "Onium_Disso_dR.h"
#include "Onium_predefine.h"
#include "predefine.h"
#include "random.h"
#include <cmath>
#include <gsl/gsl_math.h>
#include "simpleLogger.h"


// Q + Qbar --> QQbar[1S] + g
double RV1S_reco_gluon(double x[3], std::size_t dim, void * params_){
   double * params = static_cast<double *>(params_);
   double v = x[0];
   double T = x[1];
   double p = x[2];
   double M = params[0];
   double E1S = params[1];
   double a_B = params[2];
   double q = p*p/M + E1S;
   double reco = prefactor_reco_gluon * pow(q,3) * Matrix1S(p, a_B);

   if (v < 1e-4){
       double enhencement = nBplus1(q/T);
       return 2.* reco * enhencement;
   }
   else{
       double gamma = 1./std::sqrt(1.-v*v);
       double k1 = gamma*(1.+v)/T, k2 = gamma*(1.-v)/T;
       double enhencement = 2. + T/(gamma*q*v)*( fac1(q*k1) - fac1(q*k2) );
       return reco * enhencement;
   }
}

// Q + Qbar --> QQbar[2S] + g
double RV2S_reco_gluon(double x[3], std::size_t dim, void * params_){
   double * params = static_cast<double *>(params_);
   double v = x[0];
   double T = x[1];
   double p = x[2];
   double M = params[0];
   double E2S = params[1];
   double a_B = params[2];
   double q = p*p/M + E2S;
   double reco = prefactor_reco_gluon * pow(q,3) * Matrix2S(p, a_B);

   if (v < 1e-4){
       double enhencement = nBplus1(q/T);
       return 2.* reco * enhencement;
   }
   else{
       double gamma = 1./std::sqrt(1.-v*v);
       double k1 = gamma*(1.+v)/T, k2 = gamma*(1.-v)/T;
       double enhencement = 2. + T/(gamma*q*v)*( fac1(q*k1) - fac1(q*k2) );
       return reco * enhencement;
   }
}

// Q + Qbar --> QQbar[1P] + g
double RV1P_reco_gluon(double x[3], std::size_t dim, void * params_){
   double * params = static_cast<double *>(params_);
   double v = x[0];
   double T = x[1];
   double p = x[2];
   double M = params[0];
   double E1P = params[1];
   double a_B = params[2];
   double q = p*p/M + E1P;
   double reco = prefactor_reco_gluon * pow(q,3) * Matrix1P(p, a_B);

   if (v < 1e-4){
       double enhencement = nBplus1(q/T);
       return 2.* reco * enhencement;
   }
   else{
       double gamma = 1./std::sqrt(1.-v*v);
       double k1 = gamma*(1.+v)/T, k2 = gamma*(1.-v)/T;
       double enhencement = 2. + T/(gamma*q*v)*( fac1(q*k1) - fac1(q*k2) );
       return reco * enhencement;
   }
}


// Inelastic recombination, R*V
// Q + Qbar + q --> QQbar[1S] + q
double dRdp1dp2_1S_reco_ineq(double x[4], std::size_t dim, void * params_){
    double p1 = x[0];
    double c1 = x[1];
    double c2 = x[2];
    double phi = x[3];
    double * params = static_cast<double *>(params_);
    double v = params[0];
    double T = params[1];
    double p = params[2];
    double M = params[3];
    double E1S = params[4];
    double gamma = 1./std::sqrt(1.-v*v);
    double p2 = p1 + p*p/M + E1S;
    double s1 = std::sqrt(1.-c1*c1);
    double s2 = std::sqrt(1.-c2*c2);
    double phase1 = p1*nF(gamma*p1*(1.+v*c1)/T);
    double phase2 = p2*(1.-nF(gamma*p2*(1.+v*c2)/T));
    double part_angle = s1*s2*std::cos(phi)+c1*c2;
    double prop = p1*p2*(1.+part_angle)/( p1*p1+p2*p2-2.*p1*p2*part_angle );
    return phase1 * phase2 * prop; // Matrix(p, a_B) is a constant, multiply it after integration
}
// Q + Qbar + g --> QQbar[1S] + g
double dRdq1dq2_1S_reco_ineg(double x[4], std::size_t dim, void * params_){
    double q1 = x[0];
    double c1 = x[1];
    double c2 = x[2];
    double phi = x[3];
    double * params = static_cast<double *>(params_);
    double v = params[0];
    double T = params[1];
    double p = params[2];
    double M = params[3];
    double E1S = params[4];
    double gamma = 1./std::sqrt(1.-v*v);
    double q2 = q1 + p*p/M + E1S;
    double s1 = std::sqrt(1.-c1*c1);
    double s2 = std::sqrt(1.-c2*c2);
    double phase1 = q1*nB(gamma*q1*(1.+v*c1)/T);
    double phase2 = q2*(1.+nB(gamma*q2*(1.+v*c2)/T));
    double part_angle = s1*s2*std::cos(phi)+c1*c2;
    double q1q2sum = q1+q2;
    double prop = q1q2sum*q1q2sum*(1.+part_angle)/( q1*q1+q2*q2-2.*q1*q2*part_angle );
    return phase1 * phase2 * prop;  // Matrix(p, a_B) is a constant, multiply it after integration
}

// Q + Qbar + q --> QQbar[2S] + q
double dRdp1dp2_2S_reco_ineq(double x[4], std::size_t dim, void * params_){
    double p1 = x[0];
    double c1 = x[1];
    double c2 = x[2];
    double phi = x[3];
    double * params = static_cast<double *>(params_);
    double v = params[0];
    double T = params[1];
    double p = params[2];
    double M = params[3];
    double E2S = params[4];
    double gamma = 1./std::sqrt(1.-v*v);
    double p2 = p1 + p*p/M + E2S;
    double s1 = std::sqrt(1.-c1*c1);
    double s2 = std::sqrt(1.-c2*c2);
    double phase1 = p1*nF(gamma*p1*(1.+v*c1)/T);
    double phase2 = p2*(1.-nF(gamma*p2*(1.+v*c2)/T));
    double part_angle = s1*s2*std::cos(phi)+c1*c2;
    double prop = p1*p2*(1.+part_angle)/( p1*p1+p2*p2-2.*p1*p2*part_angle );
    return phase1 * phase2 * prop; // Matrix(p, a_B) is a constant, multiply it after integration
}
// Q + Qbar + g --> QQbar[2S] + g
double dRdq1dq2_2S_reco_ineg(double x[4], std::size_t dim, void * params_){
    double q1 = x[0];
    double c1 = x[1];
    double c2 = x[2];
    double phi = x[3];
    double * params = static_cast<double *>(params_);
    double v = params[0];
    double T = params[1];
    double p = params[2];
    double M = params[3];
    double E2S = params[4];
    double gamma = 1./std::sqrt(1.-v*v);
    double q2 = q1 + p*p/M + E2S;
    double s1 = std::sqrt(1.-c1*c1);
    double s2 = std::sqrt(1.-c2*c2);
    double phase1 = q1*nB(gamma*q1*(1.+v*c1)/T);
    double phase2 = q2*(1.+nB(gamma*q2*(1.+v*c2)/T));
    double part_angle = s1*s2*std::cos(phi)+c1*c2;
    double q1q2sum = q1+q2;
    double prop = q1q2sum*q1q2sum*(1.+part_angle)/( q1*q1+q2*q2-2.*q1*q2*part_angle );
    return phase1 * phase2 * prop;  // Matrix(p, a_B) is a constant, multiply it after integration
}

// Q + Qbar + q --> QQbar[1P] + q
double dRdp1dp2_1P_reco_ineq(double x[4], std::size_t dim, void * params_){
    double p1 = x[0];
    double c1 = x[1];
    double c2 = x[2];
    double phi = x[3];
    double * params = static_cast<double *>(params_);
    double v = params[0];
    double T = params[1];
    double p = params[2];
    double M = params[3];
    double E1P = params[4];
    double gamma = 1./std::sqrt(1.-v*v);
    double p2 = p1 + p*p/M + E1P;
    double s1 = std::sqrt(1.-c1*c1);
    double s2 = std::sqrt(1.-c2*c2);
    double phase1 = p1*nF(gamma*p1*(1.+v*c1)/T);
    double phase2 = p2*(1.-nF(gamma*p2*(1.+v*c2)/T));
    double part_angle = s1*s2*std::cos(phi)+c1*c2;
    double prop = p1*p2*(1.+part_angle)/( p1*p1+p2*p2-2.*p1*p2*part_angle );
    return phase1 * phase2 * prop; // Matrix(p, a_B) is a constant, multiply it after integration
}
// Q + Qbar + g --> QQbar[1P] + g
double dRdq1dq2_1P_reco_ineg(double x[4], std::size_t dim, void * params_){
    double q1 = x[0];
    double c1 = x[1];
    double c2 = x[2];
    double phi = x[3];
    double * params = static_cast<double *>(params_);
    double v = params[0];
    double T = params[1];
    double p = params[2];
    double M = params[3];
    double E1P = params[4];
    double gamma = 1./std::sqrt(1.-v*v);
    double q2 = q1 + p*p/M + E1P;
    double s1 = std::sqrt(1.-c1*c1);
    double s2 = std::sqrt(1.-c2*c2);
    double phase1 = q1*nB(gamma*q1*(1.+v*c1)/T);
    double phase2 = q2*(1.+nB(gamma*q2*(1.+v*c2)/T));
    double part_angle = s1*s2*std::cos(phi)+c1*c2;
    double q1q2sum = q1+q2;
    double prop = q1q2sum*q1q2sum*(1.+part_angle)/( q1*q1+q2*q2-2.*q1*q2*part_angle );
    return phase1 * phase2 * prop;  // Matrix(p, a_B) is a constant, multiply it after integration
}

// ---------------------------------- sampling -------------------------
// reco_real_gluon
double Sample_reco_gluon_costheta(double v, double T, double q){
    double gamma = 1./std::sqrt(1.-v*v);
    double y1 = q*gamma*(1.-v)/T;
    double max_value = nBplus1(y1);
    double x_try, y_try, result;
    do {
        x_try = Srandom::dist_costheta(Srandom::gen);
        y_try = q*gamma*(1.+ x_try*v)/T;
        result = nBplus1(y_try);
    } while (Srandom::rejection(Srandom::gen)*max_value > result);
    return x_try;
}

// reco_ineq
double f_p1_reco_important(double p1, void * params_){
    double * params = static_cast<double *>(params_);
    double v = params[0];
    double T = params[1];
    double Enl = params[2];
    double E_rel = params[3];   // E_rel = p_rel^2/M
    double gamma = 1./std::sqrt(1.-v*v);
    double k1 = gamma*(1-v)/T;  // k1 = gamma*(1-v)/T
    double k2 = gamma*(1+v)/T;  // k2 = gamma*(1+v)/T
    double p2 = p1+Enl+E_rel;
    double p1p2 = p1/p2;
    return ( fac2(k1*p1) - fac2(k2*p1) ) * p2 /(p1p2 + 1./p1p2 -2.0);
}

// importance sampling of p1 works much faster than the inverse function method, considering the overall efficiency;
// if use inverse function method, p1 sampling is easy, but the remaining integrand is a quadratic, divergent at large p1, efficiency of rejection is very low;
double Sample_reco_ineq_p1_important(double p1low, double p1up, double result_max, void * params_){
    double * params = static_cast<double *>(params_);
    double result_try, p1_try;
    int limiter = 0;
    do{
        limiter ++;
        p1_try = Srandom::rejection(Srandom::gen)*(p1up-p1low) + p1low;
        result_try = f_p1_reco_important(p1_try, params);
    } while(Srandom::rejection(Srandom::gen)*result_max > result_try && limiter < 5000);
    if (limiter >= 5000) LOG_WARNING << "Sample_reco_ineq_p1_important sampling exceed limits";
    return p1_try;
}

std::vector<double> Sample_reco_ineq(double v, double T, double p, double mass, double Enl, double maximum_f_p1){
    v = std::max(v, 1e-4);
    double gamma = 1./std::sqrt(1.-v*v);
    double p1_try, c1_try, s1_try, p2_try, c2_try, s2_try, phi_try, c_phi, s_phi, p1p2_try, part_angle, result_try;
    double p1low = 0.0;
    double p1up = 15.*T/std::sqrt(1.-v);
    double E_rel = p*p/mass;
    
    double * params_p1 = new double[4];
    params_p1[0] = v; //k1 = gamma*(1-v)/T
    params_p1[1] = T; //k2 = gamma*(1+v)/T
    params_p1[2] = Enl;
    params_p1[3] = E_rel;   // E_rel = p^2/M
    
    double * params_c1 = new double[2];
    params_c1[0] = v;
    params_c1[1] = gamma/T;

    int limiter = 0;    
    do{
        p1_try = Sample_reco_ineq_p1_important(p1low, p1up, maximum_f_p1, params_p1);
        c1_try = Sample_disso_ineq_cos1(p1_try, params_c1);
        p2_try = p1_try + E_rel + Enl;
        c2_try = Srandom::dist_costheta(Srandom::gen);
        s1_try = std::sqrt(1.-c1_try*c1_try);
        s2_try = std::sqrt(1.-c2_try*c2_try);
        phi_try = Srandom::dist_phi(Srandom::gen);
        p1p2_try = p1_try/p2_try;   // p1/p2
        c_phi = std::cos(phi_try);
        s_phi = std::sin(phi_try);
        part_angle = s1_try*s2_try*c_phi + c1_try*c2_try;
        result_try = (1.+part_angle)/(p1p2_try + 1./p1p2_try - 2.*part_angle)/2.0 * (p1p2_try + 1./p1p2_try - 2.0) * nFminus1(gamma*(1.+v*c2_try)*p2_try/T);
        limiter ++;
    } while(Srandom::rejection(Srandom::gen) >= result_try && limiter < 5000);
    if (limiter >= 5000) LOG_WARNING << "Q+Qbar+q --> QQbar + q, sampling exceed limits";
    std::vector<double> p1_final(3);
    std::vector<double> p2_final(3);
    std::vector<double> p_nl_final(3);
    p1_final = polar_to_cartisian2(p1_try, c1_try, s1_try, 1.0, 0.0);
    p2_final = polar_to_cartisian2(p2_try, c2_try, s2_try, c_phi, s_phi);
    p_nl_final = subtract_virtual_gluon(p1_final, p2_final);
    return p_nl_final;
}

// reco_ineg
double f_q1_reco_important(double q1, void * params_){
    double * params = static_cast<double *>(params_);
    double v = params[0];
    double T = params[1];
    double Enl = params[2];
    double E_rel = params[3];   // E_rel = p_rel^2/M
    double gamma = 1./std::sqrt(1.-v*v);
    double k1 = gamma*(1-v)/T;  // k1 = gamma*(1-v)/T
    double k2 = gamma*(1+v)/T;  // k2 = gamma*(1+v)/T
    double q2 = q1+Enl+E_rel;
    double q1q2 = q1/q2;
    return ( fac1(k2*q1) - fac1(k1*q1) ) * q2 * (1. + 4./(q1q2 + 1./q1q2 - 2.0) );
}

double Sample_reco_ineg_q1_important(double q1low, double q1up, double result_max, void * params_){
    double * params = static_cast<double *>(params_);
    double result_try, q1_try;
    int limiter = 0;
    do{
        limiter ++ ;
        q1_try = Srandom::rejection(Srandom::gen)*(q1up-q1low) + q1low;
        result_try = f_q1_reco_important(q1_try, params);
    } while(Srandom::rejection(Srandom::gen)*result_max > result_try && limiter < 5000);
    if (limiter >= 5000) LOG_WARNING << "Sample_reco_ineg_q1_important sampling exceed limits";
    return q1_try;
}

std::vector<double> Sample_reco_ineg(double v, double T, double p, double mass, double Enl, double maximum_f_q1){
    v = std::max(v, 1e-4);
    double gamma = 1./std::sqrt(1.-v*v);
    double q1_try, c1_try, s1_try, q2_try, c2_try, s2_try, phi_try, c_phi, s_phi, q1q2_try, part_angle, result_try;
    double q1low = 0.0;
    double q1up = 15.*T/std::sqrt(1.-v);
    double E_rel = p*p/mass;
    
    double * params_q1 = new double[4];
    params_q1[0] = v; //k1 = gamma*(1-v)/T
    params_q1[1] = T; //k2 = gamma*(1+v)/T
    params_q1[2] = Enl;
    params_q1[3] = E_rel;   // E_rel = p^2/M
    
    double * params_c1 = new double[2];
    params_c1[0] = v;
    params_c1[1] = gamma/T;
    
    int limiter = 0;
    do{
        q1_try = Sample_reco_ineg_q1_important(q1low, q1up, maximum_f_q1, params_q1);
        c1_try = Sample_disso_ineg_cos1(q1_try, params_c1);
        q2_try = q1_try + E_rel + Enl;
        c2_try = Srandom::dist_costheta(Srandom::gen);
        s1_try = std::sqrt(1.-c1_try*c1_try);
        s2_try = std::sqrt(1.-c2_try*c2_try);
        phi_try = Srandom::dist_phi(Srandom::gen);
        q1q2_try = q1_try/q2_try;   // q1/q2
        c_phi = std::cos(phi_try);
        s_phi = std::sin(phi_try);
        part_angle = s1_try*s2_try*c_phi + c1_try*c2_try;
        result_try = nBplus1(gamma*(1.+v*c2_try)*q2_try/T) * (1.+part_angle)/2.0 * (q1q2_try + 1./q1q2_try - 2.0) / (q1q2_try + 1./q1q2_try - 2.0*part_angle);
        limiter ++;
    } while(Srandom::rejection(Srandom::gen) >= result_try && limiter < 5000);
    if (limiter >= 5000) LOG_WARNING << "Q+Qbar+g --> QQbar + g, sampling exceed limits";
    std::vector<double> q1_final(3);
    std::vector<double> q2_final(3);
    std::vector<double> p_nl_final(3);
    q1_final = polar_to_cartisian2(q1_try, c1_try, s1_try, 1.0, 0.0);
    q2_final = polar_to_cartisian2(q2_try, c2_try, s2_try, c_phi, s_phi);
    p_nl_final = subtract_virtual_gluon(q1_final, q2_final);
    return p_nl_final;
}