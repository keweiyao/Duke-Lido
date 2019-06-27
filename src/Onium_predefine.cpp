#include "Onium_predefine.h"
#include <cmath>
#include "predefine.h"
#include "random.h"

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


// --------------------------- functions for sampling --------------------------
// real gluon radiation
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

/*
// inelastic scattering disso
// function used in importance sampling of p1
double f_p1_disso_important(double p1, void * params_){
    double * params = static_cast<double *>(params_);
    double k1 = params[0];        // k1 = gamma*(1-v)/T
    double k2 = params[1];        // k2 = gamma*(1+v)/T
    double p2 = p1-params[2];     // p2 = p1 - Enl
    if (p2 <= 0){
        return 0.0;
    }
    else{
        double p1p2 = p1/p2;
        return ( fac2(k1*p1) - fac2(k2*p1) ) * p1 /(p1p2 + 1./p1p2 -2.0);
    }
}

double Sample_disso_ineq_p1_important(double p1low, double p1up, double result_max, void * params_){ // result_max is an input
    double * params = static_cast<double *>(params_);
    double result_try, p1_try;
    do{
        p1_try = Srandom::rejection(Srandom::gen)*(p1up-p1low) + p1low;
        result_try = f_p1_disso_important(p1_try, params);
    } while(Srandom::rejection(Srandom::gen)*result_max > result_try);
    return p1_try;
}

double Sample_disso_ineq_cos1(double p1, void * params_){
    double * params = static_cast<double *>(params_);
    double y_try = Srandom::y_cdf(gen);
    double v = params[0];
    double B = params[1]*p1;    // B = gamma * p1/T
    double C = y_try * fac2(B*(1.+v)) + (1.-y_try) * fac2(B*(1.-v));
    return -(1. + std::log(std::exp(C) - 1.)/B )/v;
}

std::vector<double> Sample_disso_ineq(double v, double T, double mass, double Enl; double prel_up; double maximum){
    // maximum is the input for Sample_disso_ineq_p1_important
    v = std::max(v, 1e-4);
    double gamma = 1./std::sqrt(1.-v*v);
    double p1_try, c1_try, s1_try, p2_try, c2_try, s2_try, phi_try, c_phi, s_phi, p_rel, p1p2, part_angle, result_try, p1p2_try;
    double p1low = Enl;
    double p1up = 15.*T/std::sqrt(1.-v);
    
    double * params_p1 = new double[3];
    params_p1[0] = gamma*(1.-v)/T; //k1 = gamma*(1-v)/T
    params_p1[1] = gamma*(1.+v)/T; //k2 = gamma*(1+v)/T
    params_p1[2] = Enl
    
    double * params_c1 = new double[2];
    params_c1[0] = v;
    params_c1[1] = gamma/T;
    
    do{
        do{
            p_rel = Srandom::sample_inel(Srandom::gen)*prel_up;
        } while(Srandom::rejection(Srandom::gen)*max_p2Matrix1S > p2Matrix1S(p_rel));
        
        p1_try = Sample_disso_ineq_p1_important(p1low, p1up, maximum, params_p1);  // give the maximum as result_max to Sample_disso_ineq_p1_important
        p2_try = p1_try - Enl - p_rel*p_rel/mass;
        if (p2_try <= 0.0){
            result_try = 0.0;
        }
        else{
            c1_try = Sample_disso_ineq_cos1(p1_try, params_c1);
            c2_try = Srandom::dist_costheta(Srandom::gen);
            //if (v > 0.99){ c2_try = (c2_try + 1.)*4.0/gamma - 1.0; }
            s1_try = std::sqrt(1.-c1_try*c1_try);
            s2_try = std::sqrt(1.-c2_try*c2_try);
            phi_try = Srandom::dist_phi(Srandom::gen);
            p1p2_try = p1_try/p2_try;   // p1_try/p2_try
            p1p2 = (p1_try-E1S)/p1_try;
            c_phi = std::cos(phi_try);
            s_phi = std::sin(phi_try);
            part_angle = s1_try*s2_try*c_phi + c1_try*c2_try;
            result_try = (1.+part_angle)/(p1p2_try + 1./p1p2_try - 2.*part_angle)/2.0*(p1p2 + 1./p1p2 - 2.0) * nFminus1(gamma*(1.+v*c2_try)*p2_try/T)/p1p2_try;
        }
    } while(Srandom::rejection(Srandom::gen) >= result_try);
    std::vector<double> p1_final(3);
    std::vector<double> p2_final(3);
    std::vector<double> p_rel_final(3);
    std::vector<double> pQpQbar_final(6);
    double cos_rel, phi_rel;
    cos_rel = sample_cos(gen);
    phi_rel = sample_inel(gen)*TwoPi;
    p1_final = polar_to_cartisian2(p1_try, c1_try, s1_try, 1.0, 0.0);
    p2_final = polar_to_cartisian2(p2_try, c2_try, s2_try, c_phi, s_phi);
    p_rel_final = polar_to_cartisian1(p_rel, cos_rel, phi_rel);
    pQpQbar_final = add_virtual_gluon(p1_final, p2_final, p_rel_final);
    return pQpQbar_final;
}


std::vector<double> S1S_decay_ineq_test(double v, double T, double maximum){
    v = std::max(v, small_number);
    double gamma = 1./std::sqrt(1.-v*v);
    double p1_try, c1_try, s1_try, p2_try, c2_try, s2_try, phi_try, c_phi, s_phi, p_rel, p1p2, part_angle, result_try, p1p2_try;
    double p1low = E1S;
    double p1up = 15.*T/std::sqrt(1.-v);

    double * params_p1 = new double[2];
    params_p1[0] = gamma*(1.-v)/T; //k1 = gamma*(1-v)/T
    params_p1[1] = gamma*(1.+v)/T; //k2 = gamma*(1+v)/T

    double * params_c1 = new double[2];
    params_c1[0] = v;
    params_c1[1] = gamma/T;

    do{
        do{
            p_rel = sample_inel(gen)*p_1Ssam;
        } while(rejection(gen)*max_p2Matrix1S > p2Matrix1S(p_rel));

        p1_try = S1S_decay_ineq_p1_important(p1low, p1up, maximum, params_p1);
        p2_try = p1_try - E1S - p_rel*p_rel/M;
        if (p2_try <= 0.0){
            result_try = 0.0;
        }
        else{
        c1_try = S1S_decay_ineq_cos1(p1_try, params_c1);
        c2_try = sample_cos(gen);
        //if (v > 0.99){ c2_try = (c2_try + 1.)*4.0/gamma - 1.0; }
        s1_try = std::sqrt(1.-c1_try*c1_try);
        s2_try = std::sqrt(1.-c2_try*c2_try);
        phi_try = sample_inel(gen)*TwoPi;
        p1p2_try = p1_try/p2_try;   // p1_try/p2_try
        p1p2 = (p1_try-E1S)/p1_try;
        c_phi = std::cos(phi_try);
        s_phi = std::sin(phi_try);
        part_angle = s1_try*s2_try*c_phi + c1_try*c2_try;
        result_try = (1.+part_angle)/(p1p2_try + 1./p1p2_try - 2.*part_angle)/2.0*(p1p2 + 1./p1p2 - 2.0) * nFminus1(gamma*(1.+v*c2_try)*p2_try/T)/p1p2_try;
        }
    } while(rejection(gen) >= result_try);

    std::vector<double> p1p2_test(5);
    p1p2_test[0] = p1_try;
    p1p2_test[1] = c1_try;
    p1p2_test[2] = p2_try;
    p1p2_test[3] = c2_try;
    p1p2_test[4] = phi_try;
    return p1p2_test;
}
*/
