#ifndef ONIUM_DISSO_DR_H
#define ONIUM_DISSO_DR_H
#include <cmath>
#include <vector>
//-----------Quarkonium dissociation differential rate-------------------
double dRdq_1S_gluon(double q, void * params_);
double dRdq_2S_gluon(double q, void * params_);
double dRdq_1P_gluon(double q, void * params_);

double dRdp1dp2_1S_decay_ineq(double x[5], std::size_t dim, void * params_);
double dRdq1dq2_1S_decay_ineg(double x[5], std::size_t dim, void * params_);

double dRdp1dp2_2S_decay_ineq(double x[5], std::size_t dim, void * params_);
double dRdq1dq2_2S_decay_ineg(double x[5], std::size_t dim, void * params_);

double dRdp1dp2_1P_decay_ineq(double x[5], std::size_t dim, void * params_);
double dRdq1dq2_1P_decay_ineg(double x[5], std::size_t dim, void * params_);

//----------------Quarkonium dissociation sampling---------------------
// real gluon
double disso_gluon_costheta(double q, double v, double T, double random);
double Sample_reco_gluon_costheta(double v, double T, double q);

// ineq
double f_p1_disso_important(double p1, void * params_);
double Sample_disso_ineq_p1_important(double p1low, double p1up, double result_max, void * params_);
double Sample_disso_ineq_cos1(double p1, void * params_);
std::vector<double> Sample_disso_ineq(
       double v, double T, double mass, double Enl, 
       double a_B, double prel_up, double maximum_f_p1, double max_p2Matrix, 
       double(*f)(double _prel, double _aB)
);
#endif
