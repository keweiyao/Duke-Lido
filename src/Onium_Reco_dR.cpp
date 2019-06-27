#include "Onium_Reco_dR.h"
#include "Onium_predefine.h"
#include "predefine.h"
#include <cmath>
#include <gsl/gsl_math.h>

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
   double reco = prefactor_reco_gluon * pow(q,3) 
               * Matrix1S(p, a_B);

   // convert GeV^-2 to GeV * fm^3
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
   double reco = prefactor_reco_gluon * pow(q,3) 
                * Matrix2S(p, a_B);

   // convert GeV^-2 to GeV * fm^3
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
   double reco = prefactor_reco_gluon * pow(q,3) 
                 * Matrix1P(p, a_B);

   // convert GeV^-2 to GeV * fm^3
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
