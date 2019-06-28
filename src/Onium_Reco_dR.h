#ifndef ONIUM_RECO_DR_H
#define ONIUM_RECO_DR_H
#include <cmath>
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
double f_p1_reco_important(double p1, void * params_);
double Sample_reco_ineq_p1_important(double p1low, double p1up, double result_max, void * params_);
std::vector<double> S1S_reco_ineq(double v, double T, double mass, double Enl, double a_B, double prel_up, double maximum_f_p1, double max_p2Matrix, double(*f)(double _prel, double _aB);
#endif
