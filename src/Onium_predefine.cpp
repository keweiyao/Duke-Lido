#include "Onium_predefine.h"
#include <cmath>
#include "predefine.h"

const double accuracy = 0.0001;

const double li2_minus1 = -0.822467;

double nB(double x){
    double exp_factor = std::exp(-x);
    return exp_factor/(1.0 - exp_factor);
}

double nBplus1(double x){
    return 1./(1. - std::exp(-x));
}

double nF(double x){
    double exp_factor = std::exp(-x);
    return exp_factor/(1.0 + exp_factor);
}

double nFminus1(double x){
    // it is 1-n_F
    return 1./(1. + std::exp(-x));
}

double fac1(double z){
    return std::log(1. - std::exp(-z));
}

double fac2(double z){
    return std::log(1. + std::exp(-z));
}

double Li2(double z){
    int n = 0;
    double result = 0.;
    double nterm;
    do{
        n += 1;
        nterm = pow(z, n)/(n*n);
        result += nterm;
    } while(std::abs(nterm/result) > accuracy);
    return result;
}

// dipole matrix elements, Coulomb potential, free octet
const double fix_alpha_s = 0.3;
const double pot_alpha_s = 0.4;
const double alpha_s_sqd = fix_alpha_s*fix_alpha_s;

//const double MHQ1SS = 4.6500; //[GeV] bottom mass in 1S scheme

const double prefactor_Matrix1S = pow(2., 10.) * M_PI;
const double prefactor_Matrix2S = pow(2., 21.) * M_PI;
const double prefactor_Matrix1P = pow(2., 16.)/3. * M_PI;
const double prefactor_disso_gluon = 2.*fix_alpha_s/9./M_PI/M_PI; // need T/gamma^2/v
const double prefactor_disso_gluon_small_v = 1/3./M_PI/M_PI*fix_alpha_s*CF; // need 1/gamma
const double prefactor_reco_gluon = 2.*fix_alpha_s* TF / (3.*Nc); // RV prefactor

const double prefactor_disso_ineq = alpha_s_sqd*4./9./pow(M_PI,4);
// above line: for inelastic quark dissociation, spin*2, color*3, antiparticle*2, flavor*2
const double prefactor_reco_ineq = alpha_s_sqd/9./M_PI/M_PI;
// above line has 1/(Nc^2-1) * 2 * 3 * 2 * 2 = 3
const double prefactor_disso_ineg = alpha_s_sqd/6./pow(M_PI,4);
// above line: gluon inelastic decay, sum over gluon color and polarizations
const double prefactor_reco_ineg = alpha_s_sqd/24./M_PI/M_PI;
// above line: gluon inelastic reco


double get_aB(double M){
    return 2./pot_alpha_s/CF/M;
}

double get_onium_Enl_Coulomb(double M, int n, int l){
    return pot_alpha_s*CF/2./get_aB(M)/std::pow(n, 2);
}

double get_onium_Enl_Confining(double M, int n, int l, double confine_asymp_pot){
    return pot_alpha_s*CF/2./get_aB(M)/std::pow(n, 2) + confine_asymp_pot;
}

double get_onium_Mass(double MQ, int n, int l){
    return 2*MQ - get_onium_Enl_Coulomb(MQ, n, l);
}



// Matrix-elements
double Matrix1S(double prel, double a_B){
    double pa2 = pow(a_B*prel, 2.);
    return prefactor_Matrix1S * pow(a_B, 5.) * pa2/pow(1.+pa2, 6.);
}

double Matrix2S(double prel, double a_B){
    double pa2 = pow(a_B*prel, 2.);
    return prefactor_Matrix2S * pow(a_B, 5.) * pa2 * pow(1.-2.*pa2, 2.)/pow(1.+4.*pa2, 8.);
}

double Matrix1P(double prel, double a_B){
    double pa2 = pow(a_B*prel, 2.);
    return prefactor_Matrix1P * pow(a_B, 5.) * (1. - 16.*pa2 + 208.*pa2*pa2)/pow(1.+4.*pa2, 8.);
}

double p2Matrix1S(double prel, double a_B){
    return prel*prel*Matrix1S(prel, a_B);
}

double p2Matrix2S(double prel, double a_B){
    return prel*prel*Matrix2S(prel, a_B);
}

double p2Matrix1P(double prel, double a_B){
    return prel*prel*Matrix1P(prel, a_B);
}

