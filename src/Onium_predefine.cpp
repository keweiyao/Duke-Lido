#include "Onium_predefine.h"
#include <cmath>
#include "predefine.h"
#include "random.h"

const double accuracy = 0.0001;
const double li2_minus1 = -0.822467;

// find maximum of a function
// only works for positive-function with one local maximum within [xL, xH]
double find_max_noparams(double(*f)(double, double), double xL_, double xR_, double aB){
    double dfL, dfR, dfM, xM, xL = xL_, xR = xR_, dx, fM, fL, fR;
    int count = 1;
    dx = (xR-xL)/100.;
    fL = f(xL, aB);
    fR = f(xR, aB);
    dfL = f(xL+dx, aB) - fL;
    dfR = fR - f(xR-dx, aB);
    if (dfL*dfR <= 0.0){
        do{
            count += 1;
            xM = (xL+xR)/2.;
            dfM = f(xM+dx, aB) - f(xM, aB);
            fM = f(xM, aB);
            if (dfL*dfM < 0 or dfM*dfR > 0) {xR = xM; dfR = dfM;}
            else {xL = xM; dfL = dfM;}
            dx = (xR-xL)/100.;
        }while ( std::abs(dfM/dx) > accuracy and count < 10);
        return fM;
    }
    else{
        if (fL > fR) {return fL;}
        else {return fR;}
    }
}

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
const double TwoPi = 2.*M_PI;
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
double get_prel_up(double a_B, int n){
    double prel_up;
    if (n == 1){ prel_up = 4.0/a_B;}
    if (n == 2){ prel_up = 2.0/a_B;}
    return prel_up;
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

double get_max_p2Matrix(double prel_up, int n, int l, double aB){
    if (n==1){
        return find_max_noparams(&p2Matrix1S, 1e-4, prel_up, aB);
    }
    else if (n==2 && l==0){
        return find_max_noparams(&p2Matrix2S, 1e-4, prel_up, aB);
    }
    else{
        return find_max_noparams(&p2Matrix1P, 1e-4, prel_up, aB);
    }
}


// Transformations
// given three momentum and mass, calculate energy
double momentum_to_energy(double px, double py, double pz, double mass_){
    return std::sqrt(px*px+py*py+pz*pz+mass_*mass_);
}

// polar to cartesian
std::vector<double> polar_to_cartisian1(double length, double cos, double phi){
    std::vector<double> result(3);
    double sin = std::sqrt(1.-cos*cos);
    double c_phi = std::cos(phi);
    double s_phi = std::sin(phi);
    result[0] = length * sin * c_phi;
    result[1] = length * sin * s_phi;
    result[2] = length * cos;
    return result;
}

std::vector<double> polar_to_cartisian2(double length, double cos, double sin, double c_phi, double s_phi){
    std::vector<double> result(3);
    result[0] = length * sin * c_phi;
    result[1] = length * sin * s_phi;
    result[2] = length * cos;
    return result;
}

// add real gluon momentum to the decay products: QQbar
std::vector<double> add_real_gluon(std::vector<double> momentum_add, std::vector<double> momentum_rel){
    double q_x = 0.5*momentum_add[0];
    double q_y = 0.5*momentum_add[1];
    double q_z = 0.5*momentum_add[2];
    std::vector<double> pQpQbar(6);
    pQpQbar[0] = q_x + momentum_rel[0];
    pQpQbar[1] = q_y + momentum_rel[1];
    pQpQbar[2] = q_z + momentum_rel[2];
    pQpQbar[3] = q_x - momentum_rel[0];
    pQpQbar[4] = q_y - momentum_rel[1];
    pQpQbar[5] = q_z - momentum_rel[2];
    return pQpQbar;
}

// add virtual gluon momentum to the decay products: QQbar
std::vector<double> add_virtual_gluon(std::vector<double> momentum_1, std::vector<double> momentum_2, std::vector<double> momentum_rel){
    double q_x = 0.5*(momentum_1[0]-momentum_2[0]);
    double q_y = 0.5*(momentum_1[1]-momentum_2[1]);
    double q_z = 0.5*(momentum_1[2]-momentum_2[2]);
    std::vector<double> pQpQbar(6);
    pQpQbar[0] = q_x + momentum_rel[0];
    pQpQbar[1] = q_y + momentum_rel[1];
    pQpQbar[2] = q_z + momentum_rel[2];
    pQpQbar[3] = q_x - momentum_rel[0];
    pQpQbar[4] = q_y - momentum_rel[1];
    pQpQbar[5] = q_z - momentum_rel[2];
    return pQpQbar;
}

// subtract real gluon momentum to the reco produce: |nlm>
std::vector<double> subtract_real_gluon(std::vector<double> momentum_subtract){
    std::vector<double> p_nl(3);    // nl here represents quarkonia
    p_nl[0] = -momentum_subtract[0];
    p_nl[1] = -momentum_subtract[1];
    p_nl[2] = -momentum_subtract[2];
    return p_nl;
}

// subtract virtual gluon momentum to the reco produce: |nlm>
std::vector<double> subtract_virtual_gluon(std::vector<double> momentum_1, std::vector<double> momentum_2){
    std::vector<double> p_nl(3);
    p_nl[0] = momentum_1[0]-momentum_2[0];
    p_nl[1] = momentum_1[1]-momentum_2[1];
    p_nl[2] = momentum_1[2]-momentum_2[2];
    return p_nl;
}


