#ifndef ONIUM_RECO_DR_H
#define ONIUM_RECO_DR_H
#include <cmath>
#include <vector>
// Q + Qbar --> QQbar[1S] + g
double RV1S_reco_gluon(double x[3], std::size_t dim, void * params_);
// Q + Qbar --> QQbar[2S] + g
double RV2S_reco_gluon(double x[3], std::size_t dim, void * params_);

// Q + Qbar --> QQbar[1P] + g
double RV1P_reco_gluon(double x[3], std::size_t dim, void * params_);


// Inelastic recombination, R*V
// Q + Qbar + q --> QQbar[1S] + q
double dRdp1dp2_1S_reco_ineq(double x[4], std::size_t dim, void * params_);
// Q + Qbar + g --> QQbar[1S] + g
double dRdq1dq2_1S_reco_ineg(double x[4], std::size_t dim, void * params_);
// Q + Qbar + q --> QQbar[2S] + q
double dRdp1dp2_2S_reco_ineq(double x[4], std::size_t dim, void * params_);
// Q + Qbar + g --> QQbar[2S] + g
double dRdq1dq2_2S_reco_ineg(double x[4], std::size_t dim, void * params_);
// Q + Qbar + q --> QQbar[1P] + q
double dRdp1dp2_1P_reco_ineq(double x[4], std::size_t dim, void * params_);
// Q + Qbar + g --> QQbar[1P] + g
double dRdq1dq2_1P_reco_ineg(double x[4], std::size_t dim, void * params_);


// Sampling
// real_gluon
double Sample_reco_gluon_costheta(double v, double T, double q);

// ineq
double f_p1_reco_important(double p1, void * params_);
double Sample_reco_ineq_p1_important(double p1low, double p1up, double result_max, void * params_);
std::vector<double> Sample_reco_ineq(
           double v, double T, double p, 
           double mass, double Enl, double maximum_f_p1
                                     );

//ineg
double f_q1_reco_important(double q1, void * params_);
double Sample_reco_ineg_q1_important(double q1low, double q1up, double result_max, void * params_);
std::vector<double> Sample_reco_ineg(
            double v, double T, double p,
            double mass, double Enl, double maximum_f_q1
                                     );
#endif
