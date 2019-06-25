#ifndef ONIUM_DISSO_DR_H
#define ONIUM_DISSO_DR_H
#include <cmath>
//-----------Quarkonium dissocation differential rate-------------------
double dRdq_1S_gluon(double q, void * params_);
double dRdq_2S_gluon(double q, void * params_);
double dRdq_1P_gluon(double q, void * params_);

double dRdp1dp2_1S_decay_ineq(double x[5], std::size_t dim, void * params_);
double dRdq1dq2_1S_decay_ineg(double x[5], std::size_t dim, void * params_);

double dRdp1dp2_2S_decay_ineq(double x[5], std::size_t dim, void * params_);
double dRdq1dq2_2S_decay_ineg(double x[5], std::size_t dim, void * params_);

double dRdp1dp2_1P_decay_ineq(double x[5], std::size_t dim, void * params_);
double dRdq1dq2_1P_decay_ineg(double x[5], std::size_t dim, void * params_);
#endif
