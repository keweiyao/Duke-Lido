#include "Onium_Rate.h"
#include "integrator.h"
#include "sampler.h"
#include "minimizer.h"
#include "random.h"
#include "approx_functions.h"
#include "matrix_elements.h"
#include "Langevin.h"
#include "Onium_Disso_dR.h"
#include "Onium_Reco_dR.h"
#include "Onium_predefine.h"
#include <iostream>
#include <cmath>
#include <gsl/gsl_math.h>


// --------------------------- Quarkonium 2->2 disso rate: QQbar[nl] + g --> Q + Qbar
template<>
OniumDissoRate22<2, double(*)(double, void *)>::
    OniumDissoRate22(std::string Name, std::string configfile, int n, int l, double(*f)(double, void *)):
StochasticBase<2>(Name+"/rate", configfile),
_f(f)
{
   // read configfile
   boost::property_tree::ptree config;
   std::ifstream input(configfile);
   read_xml(input, config);

   std::vector<std::string> strs;
   boost::split(strs, Name, boost::is_any_of("/"));
   std::string model_name = strs[0];
   std::string process_name = strs[1];
   auto tree = config.get_child(model_name + "." + process_name);
   _mass = tree.get<double>("mass");
   _n = n;
   _l = l;
   _Enl = get_onium_Enl_Coulomb(_mass, _n, _l);
   _aB = get_aB(_mass);
   _active = (tree.get<std::string>("<xmlattr>.status")=="active")?true:false;
}

// Calculate R = int_E[nl]^inf dR/dq * dq
template <>
scalar OniumDissoRate22<2, double(*)(double, void*)>::
calculate_scalar(std::vector<double> parameters){
    double v = parameters[0];
    double T = parameters[1];
    
    // dR/dq to be intergated
    auto dRdq = [v, T, this](double q){
        double params[5] = {v, T, this->_mass, this->_Enl, this->_aB};
        double result = this->_f(q, params);
        return result;
    };
    double qmin = _Enl;
    double qmax = 100*_Enl;
    double error;
    double res = quad_1d(dRdq, {qmin, qmax}, error);
    return scalar{res}; // prefactor_disso_gluon has been accounted for in differential rate
}
// Find the max of dR/dq
template <>
scalar OniumDissoRate22<2, double(*)(double, void*)>::
find_max(std::vector<double> parameters){
    double v = parameters[0];
    double T = parameters[1];
    auto dRdq = [v, T, this](double q){
        double params[5] = {v, T, this->_mass, this->_Enl, this->_aB};
        double result = this->_f(q, params);
        return -result;
    };
    double qmin = _Enl;
    double qmax = 100*_Enl;
    double res = -minimize_1d(dRdq, {qmin, qmax}, 1e-3, 100, 100);
    return scalar{res*1.2};
}

// Sample Final states using dR/dq
template<>
void OniumDissoRate22<2, double(*)(double, void *)>::
sample(std::vector<double> parameters, std::vector< fourvec > & FS){
    double v = parameters[0];
    double T = parameters[1];
    auto dRdq = [v, T, this](double q){
        double params[5] = {v, T, this->_mass, this->_Enl, this->_aB};
        double result = this->_f(q, params);
        return result;
    };
    double qmin = _Enl;
    double qmax = 100*_Enl;
    // sample gluon energy
    double q = sample_1d(dRdq, {qmin, qmax},StochasticBase<2>::GetFmax(parameters).s);
    // sample relative momentum of Q and Qbar
    double p_rel = std::sqrt((q-_Enl)*_mass);
    // sample costheta of q
    double r = Srandom::rejection(Srandom::gen);
    double cos = disso_gluon_costheta(q, v, T, r);
    // sample phi of q
    double phi = Srandom::dist_phi(Srandom::gen);
    // sample costheta and phi of p_rel
    double cos_rel = Srandom::dist_costheta(Srandom::gen);
    double phi_rel = Srandom::dist_phi(Srandom::gen);
    
    std::vector<double> momentum_gluon(3);
    std::vector<double> momentum_rel(3);
    std::vector<double> pQpQbar_final(6);
    momentum_gluon = polar_to_cartisian1(q, cos, phi);
    momentum_rel = polar_to_cartisian1(p_rel, cos_rel, phi_rel);
    pQpQbar_final = add_real_gluon(momentum_gluon, momentum_rel);
    double E_Q = momentum_to_energy(_mass, pQpQbar_final[0], pQpQbar_final[1], pQpQbar_final[2]);
    double E_Qbar = momentum_to_energy(_mass, pQpQbar_final[3], pQpQbar_final[4], pQpQbar_final[5]);
    FS.clear();
    FS.resize(2);
    // In Rest frame
    FS[0] = fourvec{E_Q, pQpQbar_final[0], pQpQbar_final[1], pQpQbar_final[2]};
    FS[1] = fourvec{E_Qbar, pQpQbar_final[3], pQpQbar_final[4], pQpQbar_final[5]};
    // boost to hydrocell frame
    FS[0] = FS[0].boost_back(0,0,v);
    FS[1] = FS[1].boost_back(0,0,v);
}


template <size_t N, typename F>
fourvec OniumDissoRate22<N, F>::calculate_fourvec(std::vector<double> parameters){
	return fourvec::unity();
}
template <size_t N, typename F>
tensor OniumDissoRate22<N, F>::calculate_tensor(std::vector<double> parameters){
	return tensor::unity();
}


// ------------- Quarkonium 2->2 reco rate: Q + Qbar --> QQbar[nl] + g
template<>
OniumRecoRate22<3, double(*)(double *, std::size_t, void *)>::
OniumRecoRate22(std::string Name, std::string configfile, int n, int l, double(*f)(double *, std::size_t, void *)):
StochasticBase<3>(Name+"/rate", configfile),
_f(f)
{
    // read configfile
    boost::property_tree::ptree config;
    std::ifstream input(configfile);
    read_xml(input, config);
    
    std::vector<std::string> strs;
    boost::split(strs, Name, boost::is_any_of("/"));
    std::string model_name = strs[0];
    std::string process_name = strs[1];
    auto tree = config.get_child(model_name + "." + process_name);
    _mass = tree.get<double>("mass");
    _n = n;
    _l = l;
    _Enl = get_onium_Enl_Coulomb(_mass, _n, _l);
    _aB = get_aB(_mass);
    _active = (tree.get<std::string>("<xmlattr>.status")=="active")?true:false;
}

// 2->2 reco rate
template <>
scalar OniumRecoRate22<3, double(*)(double *, std::size_t, void *)>::
find_max(std::vector<double> parameters){
    return scalar::unity();
}


// 2->2 reco rate
template <>
scalar OniumRecoRate22<3, double(*)(double *, std::size_t, void *)>::
calculate_scalar(std::vector<double> parameters){
    double x[3];
    x[0] = parameters[0]; // v
    x[1] = parameters[1]; // T
    x[2] = std::exp(parameters[2]); // exp(ln_p_rel)
    double params[3] = {this->_mass, this->_Enl, this->_aB};
    double result = prefactor_reco_gluon * this->_f(x, 3, params);
    return scalar{result};
}


// sample reco_real_gluon
template<>
void OniumRecoRate22<3, double(*)(double *, std::size_t, void *)>::
sample(std::vector<double> parameters, std::vector< fourvec > & FS){
    double v = parameters[0];
    double T = parameters[1];
    double p = std::exp(parameters[2]);      // relative momentum in QQbar rest frame
    
    double q = p*p/_mass + _Enl; // emitted gluon energy
    double cos = Sample_reco_gluon_costheta(v, T, q); // emitted gluon direction
    double phi = Srandom::dist_phi(Srandom::gen);
    
    std::vector<double> momentum_gluon(3);
    std::vector<double> momentum_onium(3);
    momentum_gluon = polar_to_cartisian1(q, cos, phi);
    momentum_onium = subtract_real_gluon(momentum_gluon);
    double E_onium = momentum_to_energy(2.*_mass-_Enl, momentum_onium[0], momentum_onium[1], momentum_onium[2]);
    FS.clear();
    FS.resize(1);
    // In Rest frame of QQbar
    FS[0] = fourvec{E_onium, momentum_onium[0], momentum_onium[1], momentum_onium[2]};
    // boost to hydrocell frame
    FS[0] = FS[0].boost_back(0,0,v);
}

template <size_t N, typename F>
fourvec OniumRecoRate22<N, F>::calculate_fourvec(std::vector<double> parameters){
    return fourvec::unity();
}
template <size_t N, typename F>
tensor OniumRecoRate22<N, F>::calculate_tensor(std::vector<double> parameters){
    return tensor::unity();
}


// ------------------------------- Quarkonium 2->3 ineq disso rate: QQbar[nl] + q --> Q + Qbar + q
template<>
OniumDissoRate23q<2, double(*)(double *, std::size_t, void *), double(*)(double, void*), double(*)(double, double)>::
OniumDissoRate23q(std::string Name, std::string configfile, int n, int l, double(*f1)(double *, std::size_t, void *), double(*f2)(double, void*), double(*f3)(double, double)):
StochasticBase<2>(Name+"/rate", configfile),
_f1(f1),
_f2(f2),
_f3(f3)
{
    // read configfile
    boost::property_tree::ptree config;
    std::ifstream input(configfile);
    read_xml(input, config);
    
    std::vector<std::string> strs;
    boost::split(strs, Name, boost::is_any_of("/"));
    std::string model_name = strs[0];
    std::string process_name = strs[1];
    auto tree = config.get_child(model_name + "." + process_name);
    _mass = tree.get<double>("mass");
    _n = n;
    _l = l;
    _Enl = get_onium_Enl_Coulomb(_mass, _n, _l);
    _aB = get_aB(_mass);
    _prel_up = get_prel_up(_aB, _n);
    _max_p2Matrix = get_max_p2Matrix(_prel_up, _n, _l, _aB);
    _active = (tree.get<std::string>("<xmlattr>.status")=="active")?true:false;
}

// Calculate R = int dp1 dp2 dR/dp1/dp2
template <>
scalar OniumDissoRate23q<2, double(*)(double *, std::size_t, void *), double(*)(double, void*), double(*)(double, double)>::
calculate_scalar(std::vector<double> parameters){
    double v = parameters[0];
    double T = parameters[1];
    
    // dR/dp1/dp2 to be intergated
    auto dRdp1dp2 = [v, T, this](double * x){
        double params[5] = {v, T, this->_mass, this->_Enl, this->_aB};
        double result = this->_f1(x, 5, params);
        return result;
    };
    double p1up = 15.*T/std::sqrt(1.-v);

    double xl[5] = { _Enl, -1., 0., -1., 0. };
    double xu[5] = { p1up, 1., _prel_up, 1., TwoPi };
    double error;
    double res = prefactor_disso_ineq * vegas(dRdp1dp2, 5, xl, xu, error, 5000);
    return scalar{res};
}

// find max of f_p1
template <>
scalar OniumDissoRate23q<2, double(*)(double*, std::size_t, void*), double(*)(double, void*), double(*)(double, double)>::
find_max(std::vector<double> parameters){
    double v = parameters[0];
    double T = parameters[1];
    auto f_p1 = [v, T, this](double p1){
        double params[3] = {v, T, this->_Enl};
        double result = this->_f2(p1, params);
        return -result;
    };
    double p1_min = 0.0;
    double p1_max = 15.*T/std::sqrt(1.-v);
    double res = -minimize_1d(f_p1, {p1_min, p1_max}, 1e-3, 100, 100);
    return scalar{res*1.2};
}

// Sample final states using dR/dp1/dp2
template<>
void OniumDissoRate23q<2, double(*)(double *, std::size_t, void *), double(*)(double, void*), double(*)(double, double)>::
sample(std::vector<double> parameters, std::vector< fourvec > & FS){
    double v = parameters[0];
    double T = parameters[1];
    std::vector<double> pQpQbar_final = Sample_disso_ineq(v, T, _mass, _Enl, _aB, _prel_up, StochasticBase<2>::GetFmax(parameters).s, _max_p2Matrix, this->_f3); //vector function: v, T, mQ, Enl, aB, prel_up, I_p1_max, max_p2Matrix, p2Matrix(prel, aB)

    double E_Q = momentum_to_energy(_mass, pQpQbar_final[0], pQpQbar_final[1], pQpQbar_final[2]);
    double E_Qbar = momentum_to_energy(_mass, pQpQbar_final[3], pQpQbar_final[4], pQpQbar_final[5]);
    FS.clear();
    FS.resize(2);
    // In Rest frame
    FS[0] = fourvec{E_Q, pQpQbar_final[0], pQpQbar_final[1], pQpQbar_final[2]};
    FS[1] = fourvec{E_Qbar, pQpQbar_final[3], pQpQbar_final[4], pQpQbar_final[5]};
    // boost to hydrocell frame
    FS[0] = FS[0].boost_back(0,0,v);
    FS[1] = FS[1].boost_back(0,0,v);
}

template <size_t N, typename F1, typename F2, typename F3>
fourvec OniumDissoRate23q<N, F1, F2, F3>::calculate_fourvec(std::vector<double> parameters){
    return fourvec::unity();
}
template <size_t N, typename F1, typename F2, typename F3>
tensor OniumDissoRate23q<N, F1, F2, F3>::calculate_tensor(std::vector<double> parameters){
    return tensor::unity();
}

// ------------------------------- Quarkonium 3->2 ineq reco rate: Q + Qbar + q --> QQbar[nl] + q
template<>
OniumRecoRate32q<3, double(*)(double *, std::size_t, void *), 
                    double(*)(double, void *), 
                    double(*)(double, double)>::
OniumRecoRate32q(std::string Name, std::string configfile, 
                 int n, int l, 
                 double(*f1)(double *, std::size_t, void *), 
                 double(*f2)(double, void *), 
                 double(*f3)(double, double)):
StochasticBase<3>(Name+"/rate", configfile),
_f1(f1),
_f2(f2),
_f3(f3)
{
    // read configfile
    boost::property_tree::ptree config;
    std::ifstream input(configfile);
    read_xml(input, config);
    
    std::vector<std::string> strs;
    boost::split(strs, Name, boost::is_any_of("/"));
    std::string model_name = strs[0];
    std::string process_name = strs[1];
    auto tree = config.get_child(model_name + "." + process_name);
    _mass = tree.get<double>("mass");
    _n = n;
    _l = l;
    _Enl = get_onium_Enl_Coulomb(_mass, _n, _l);
    _aB = get_aB(_mass);
    _active = (tree.get<std::string>("<xmlattr>.status")=="active")?true:false;
}

// Calculate R = int dp1 dp2 dR/dp1/dp2
template <>
scalar OniumRecoRate32q<3, 
  double(*)(double *, std::size_t, void *), 
  double(*)(double, void *), 
  double(*)(double, double)>::
calculate_scalar(std::vector<double> parameters){
    double v = parameters[0];
    double T = parameters[1];
    double p = std::exp(parameters[2]);      // relative momentum in QQbar rest frame
    
    // dR/dp1/dp2 to be intergated
    auto dRdp1dp2 = [v, T, p, this](double * x){
        double params[5] = {v, T, p, this->_mass, this->_Enl};
        double result = this->_f1(x, 4, params);
        return result;
    };
    double p1up = 15.*T/std::sqrt(1.-v);
    
    double xl[4] = { 0., -1., -1., 0. };
    double xu[4] = { p1up, 1., 1., TwoPi };
    double error;
    double res = prefactor_reco_ineq * vegas(dRdp1dp2, 4, xl, xu, error, 5000);
    res = res * this->_f3(p, _aB);
    return scalar{res};
}

// find max of f_p1
template <>
scalar OniumRecoRate32q<3, 
  double(*)(double *, std::size_t, void *), 
  double(*)(double, void *), 
  double(*)(double, double)>::
find_max(std::vector<double> parameters){
    double v = parameters[0];
    double T = parameters[1];
    double p = std::exp(parameters[2]);// relative momentum in QQbar rest frame
    auto f_p1 = [v, T, p, this](double p1){
        double params[4] = {v, T, this->_Enl, p*p/this->_mass};
        double result = this->_f2(p1, params);
        return -result;
    };
    double p1_min = 0.0;
    double p1_max = 15.*T/std::sqrt(1.-v);
    double res = -minimize_1d(f_p1, {p1_min, p1_max}, 1e-3, 100, 100);
    return scalar{res*1.2};
}

// Sample final states using dR/dp1/dp2
template<>
void OniumRecoRate32q<3, double(*)(double *, std::size_t, void *), 
                         double(*)(double, void *), 
                         double(*)(double, double)>::
sample(std::vector<double> parameters, std::vector< fourvec > & FS){
    double v = parameters[0];
    double T = parameters[1];
    double p = std::exp(parameters[2]); // relative momentum in QQbar rest frame
    std::vector<double> ponium_final = Sample_reco_ineq(v, T, p, _mass, _Enl, StochasticBase<3>::GetFmax(parameters).s);

    double E_onium = momentum_to_energy(2*_mass-_Enl, ponium_final[0], ponium_final[1], ponium_final[2]);
    FS.clear();
    FS.resize(1);
    // In Rest frame
    FS[0] = fourvec{E_onium, ponium_final[0], ponium_final[1], ponium_final[2]};
    // boost to hydrocell frame
    FS[0] = FS[0].boost_back(0,0,v);
}

template <size_t N, typename F1, typename F2, typename F3>
fourvec OniumRecoRate32q<N, F1, F2, F3>::calculate_fourvec(std::vector<double> parameters){
    return fourvec::unity();
}
template <size_t N, typename F1, typename F2, typename F3>
tensor OniumRecoRate32q<N, F1, F2, F3>::calculate_tensor(std::vector<double> parameters){
    return tensor::unity();
}


// ------------------------------- Quarkonium 2->3 ineg disso rate: QQbar[nl] + g --> Q + Qbar + g
template<>
OniumDissoRate23g<2, double(*)(double *, std::size_t, void *), double(*)(double, void*), double(*)(double, double)>::
OniumDissoRate23g(std::string Name, std::string configfile, int n, int l, double(*f1)(double *, std::size_t, void *), double(*f2)(double, void*), double(*f3)(double, double)):
StochasticBase<2>(Name+"/rate", configfile),
_f1(f1),
_f2(f2),
_f3(f3)
{
    // read configfile
    boost::property_tree::ptree config;
    std::ifstream input(configfile);
    read_xml(input, config);
    
    std::vector<std::string> strs;
    boost::split(strs, Name, boost::is_any_of("/"));
    std::string model_name = strs[0];
    std::string process_name = strs[1];
    auto tree = config.get_child(model_name + "." + process_name);
    _mass = tree.get<double>("mass");
    _n = n;
    _l = l;
    _Enl = get_onium_Enl_Coulomb(_mass, _n, _l);
    _aB = get_aB(_mass);
    _prel_up = get_prel_up(_aB, _n);
    _max_p2Matrix = get_max_p2Matrix(_prel_up, _n, _l, _aB);
    _active = (tree.get<std::string>("<xmlattr>.status")=="active")?true:false;
}

// Calculate R = int dq1 dq2 dR/dq1/dq2
template <>
scalar OniumDissoRate23g<2, double(*)(double *, std::size_t, void *), double(*)(double, void*), double(*)(double, double)>::
calculate_scalar(std::vector<double> parameters){
    double v = parameters[0];
    double T = parameters[1];
    
    // dR/dq1/dq2 to be intergated
    auto dRdq1dq2 = [v, T, this](double * x){
        double params[5] = {v, T, this->_mass, this->_Enl, this->_aB};
        double result = this->_f1(x, 5, params);
        return result;
    };
    double q1up = 15.*T/std::sqrt(1.-v);
    
    double xl[5] = { _Enl, -1., 0., -1., 0. };
    double xu[5] = { q1up, 1., _prel_up, 1., TwoPi };
    double error;
    double res = prefactor_disso_ineg * vegas(dRdq1dq2, 5, xl, xu, error, 5000);
    return scalar{res};
}

// find max of f_q1
template <>
scalar OniumDissoRate23g<2, double(*)(double*, std::size_t, void*), double(*)(double, void*), double(*)(double, double)>::
find_max(std::vector<double> parameters){
    double v = parameters[0];
    double T = parameters[1];
    auto f_q1 = [v, T, this](double q1){
        double params[3] = {v, T, this->_Enl};
        double result = this->_f2(q1, params);
        return -result;
    };
    double q1_min = 0.0;
    double q1_max = 15.*T/std::sqrt(1.-v);
    double res = -minimize_1d(f_q1, {q1_min, q1_max}, 1e-3, 100, 100);
    return scalar{res*1.2};
}

// Sample final states using dR/dq1/dq2
template<>
void OniumDissoRate23g<2, double(*)(double *, std::size_t, void *), double(*)(double, void*), double(*)(double, double)>::
sample(std::vector<double> parameters, std::vector< fourvec > & FS){
    double v = parameters[0];
    double T = parameters[1];
    std::vector<double> pQpQbar_final = Sample_disso_ineg(v, T, _mass, _Enl, _aB, _prel_up, StochasticBase<2>::GetFmax(parameters).s, _max_p2Matrix, this->_f3); //vector function: v, T, mQ, Enl, aB, prel_up, I_q1_max, max_p2Matrix, p2Matrix(prel, aB)
    
    double E_Q = momentum_to_energy(_mass, pQpQbar_final[0], pQpQbar_final[1], pQpQbar_final[2]);
    double E_Qbar = momentum_to_energy(_mass, pQpQbar_final[3], pQpQbar_final[4], pQpQbar_final[5]);
    FS.clear();
    FS.resize(2);
    // In Rest frame
    FS[0] = fourvec{E_Q, pQpQbar_final[0], pQpQbar_final[1], pQpQbar_final[2]};
    FS[1] = fourvec{E_Qbar, pQpQbar_final[3], pQpQbar_final[4], pQpQbar_final[5]};
    // boost to hydrocell frame
    FS[0] = FS[0].boost_back(0,0,v);
    FS[1] = FS[1].boost_back(0,0,v);
}

template <size_t N, typename F1, typename F2, typename F3>
fourvec OniumDissoRate23g<N, F1, F2, F3>::calculate_fourvec(std::vector<double> parameters){
    return fourvec::unity();
}
template <size_t N, typename F1, typename F2, typename F3>
tensor OniumDissoRate23g<N, F1, F2, F3>::calculate_tensor(std::vector<double> parameters){
    return tensor::unity();
}


// ------------------------------- Quarkonium 3->2 ineg reco rate: Q + Qbar + g --> QQbar[nl] + g
template<>
OniumRecoRate32g<3, double(*)(double *, std::size_t, void *),
double(*)(double, void *),
double(*)(double, double)>::
OniumRecoRate32g(std::string Name, std::string configfile,
                 int n, int l,
                 double(*f1)(double *, std::size_t, void *),
                 double(*f2)(double, void *),
                 double(*f3)(double, double)):
StochasticBase<3>(Name+"/rate", configfile),
_f1(f1),
_f2(f2),
_f3(f3)
{
    // read configfile
    boost::property_tree::ptree config;
    std::ifstream input(configfile);
    read_xml(input, config);
    
    std::vector<std::string> strs;
    boost::split(strs, Name, boost::is_any_of("/"));
    std::string model_name = strs[0];
    std::string process_name = strs[1];
    auto tree = config.get_child(model_name + "." + process_name);
    _mass = tree.get<double>("mass");
    _n = n;
    _l = l;
    _Enl = get_onium_Enl_Coulomb(_mass, _n, _l);
    _aB = get_aB(_mass);
    _active = (tree.get<std::string>("<xmlattr>.status")=="active")?true:false;
}

// Calculate R = int dq1 dq2 dR/dq1/dq2
template <>
scalar OniumRecoRate32g<3,
double(*)(double *, std::size_t, void *),
double(*)(double, void *),
double(*)(double, double)>::
calculate_scalar(std::vector<double> parameters){
    double v = parameters[0];
    double T = parameters[1];
    double p = std::exp(parameters[2]);      // relative momentum in QQbar rest frame
    
    // dR/dq1/dq2 to be intergated
    auto dRdq1dq2 = [v, T, p, this](double * x){
        double params[5] = {v, T, p, this->_mass, this->_Enl};
        double result = this->_f1(x, 4, params);
        return result;
    };
    double q1up = 15.*T/std::sqrt(1.-v);
    
    double xl[4] = { 0., -1., -1., 0. };
    double xu[4] = { q1up, 1., 1., TwoPi };
    double error;
    double res = prefactor_reco_ineg * vegas(dRdq1dq2, 4, xl, xu, error, 5000);
    res = res * this->_f3(p, _aB);
    return scalar{res};
}

// find max of f_q1
template <>
scalar OniumRecoRate32g<3,
double(*)(double *, std::size_t, void *),
double(*)(double, void *),
double(*)(double, double)>::
find_max(std::vector<double> parameters){
    double v = parameters[0];
    double T = parameters[1];
    double p = std::exp(parameters[2]);// relative momentum in QQbar rest frame
    auto f_q1 = [v, T, p, this](double q1){
        double params[4] = {v, T, this->_Enl, p*p/this->_mass};
        double result = this->_f2(q1, params);
        return -result;
    };
    double q1_min = 0.0;
    double q1_max = 15.*T/std::sqrt(1.-v);
    double res = -minimize_1d(f_q1, {q1_min, q1_max}, 1e-3, 100, 100);
    return scalar{res*1.2};
}

// Sample final states using dR/dq1/dq2
template<>
void OniumRecoRate32g<3, double(*)(double *, std::size_t, void *),
double(*)(double, void *),
double(*)(double, double)>::
sample(std::vector<double> parameters, std::vector< fourvec > & FS){
    double v = parameters[0];
    double T = parameters[1];
    double p = std::exp(parameters[2]); // relative momentum in QQbar rest frame
    std::vector<double> ponium_final = Sample_reco_ineg(v, T, p, _mass, _Enl, StochasticBase<3>::GetFmax(parameters).s);
    
    double E_onium = momentum_to_energy(2*_mass-_Enl, ponium_final[0], ponium_final[1], ponium_final[2]);
    FS.clear();
    FS.resize(1);
    // In Rest frame
    FS[0] = fourvec{E_onium, ponium_final[0], ponium_final[1], ponium_final[2]};
    // boost to hydrocell frame
    FS[0] = FS[0].boost_back(0,0,v);
}

template <size_t N, typename F1, typename F2, typename F3>
fourvec OniumRecoRate32g<N, F1, F2, F3>::calculate_fourvec(std::vector<double> parameters){
    return fourvec::unity();
}
template <size_t N, typename F1, typename F2, typename F3>
tensor OniumRecoRate32g<N, F1, F2, F3>::calculate_tensor(std::vector<double> parameters){
    return tensor::unity();
}



template class OniumDissoRate22<2,double(*)(double, void*)>;
template class OniumRecoRate22<3, double(*)(double *, std::size_t, void*)>;
template class OniumDissoRate23q<2, 
 double(*)(double *, std::size_t, void *), 
 double(*)(double, void *), 
 double(*)(double, double)>;
template class OniumRecoRate32q<3, 
 double(*)(double *, std::size_t, void *), 
 double(*)(double, void *), 
 double(*)(double, double)
>;
template class OniumDissoRate23g<2,
double(*)(double *, std::size_t, void *),
double(*)(double, void *),
double(*)(double, double)>;
template class OniumRecoRate32g<3,
double(*)(double *, std::size_t, void *),
double(*)(double, void *),
double(*)(double, double)
>;

