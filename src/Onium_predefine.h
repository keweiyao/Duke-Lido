#ifndef ONIUM_PREDEFINE_H
#define ONIUM_PREDEFINE_H

#include <vector>

extern const double accuracy;
extern const double li2_minus1;
extern const double fix_alpha_s;
extern const double pot_alpha_s;
extern const double alpha_s_sqd;

extern const double prefactor_Matrix1S;
extern const double prefactor_Matrix2S;
extern const double prefactor_Matrix1P;
extern const double prefactor_disso_gluon;
extern const double prefactor_disso_gluon_small_v;
extern const double prefactor_reco_gluon;

extern const double prefactor_disso_ineq;
extern const double prefactor_reco_ineq;
extern const double prefactor_disso_ineg;
extern const double prefactor_reco_ineg;


double nB(double x);
double nBplus1(double x);
double nF(double x);
double nFminus1(double x);
double fac1(double z);
double fac2(double z);
double Li2(double z);

// Matrix-elements
double Matrix1S(double prel, double a_B);
double Matrix2S(double prel, double a_B);
double Matrix1P(double prel, double a_B);
double p2Matrix1S(double prel, double a_B);
double p2Matrix2S(double prel, double a_B);
double p2Matrix1P(double prel, double a_B);

double get_aB(double M);
double get_onium_Enl_Coulomb(double M, int n, int l);
double get_onium_Enl_Confining(double M, int n, int l, double confine_asymp_pot);
double get_onium_Mass(double MQ, int n, int l);

// Transformations
double momentum_to_energy(double px, double py, double pz, double mass_);
std::vector<double> polar_to_cartisian1(double length, double cos, double phi);
std::vector<double> polar_to_cartisian2(double length, double cos, double sin, double c_phi, double s_phi);
std::vector<double> add_real_gluon(std::vector<double> momentum_add, std::vector<double> momentum_rel);
std::vector<double> add_virtual_gluon(std::vector<double> momentum_1, std::vector<double> momentum_2, std::vector<double> momentum_rel);
std::vector<double> subtract_real_gluon(std::vector<double> momentum_subtract);
std::vector<double> subtract_virtual_gluon(std::vector<double> momentum_1, std::vector<double> momentum_2);

// sampling
double Sample_reco_gluon_costheta(double v, double T, double q);
#endif
