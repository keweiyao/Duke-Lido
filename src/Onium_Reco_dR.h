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

#endif
