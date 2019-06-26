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


// Quarkonium 2->2 disso rate: QQbar[nl] + g --> Q + Qbar
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
	return scalar{res*1.5};
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
	return scalar{res};
}


template <size_t N, typename F>
fourvec OniumDissoRate22<N, F>::calculate_fourvec(std::vector<double> parameters){
	return fourvec::unity();
}
template <size_t N, typename F>
tensor OniumDissoRate22<N, F>::calculate_tensor(std::vector<double> parameters){
	return tensor::unity();
}

/*
// Quarkonium 2->2 reco rate: Q + Qbar --> QQbar[nl] + g
template<>
OniumRecoRate22<3, double(*)(double *, std::size_t, void *)>::
OniumRecoRate22(std::string Name, std::string configfile, int n, int l, double(*f)(double, std::size_t, void *)):
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
scalar OniumRecoRate22<3, double(*)(double, std::size_t, void*)>::
calculate_scalar(std::vector<double> parameters){
    //double v = parameters[0];
    //double T = parameters[1];
    //double p = parameters[2];
    double params[3] = {this->_mass, this->_Enl, this->_aB};
    double result = this->_f(parameters, 3, params);
    return scalar{result};
}

template <size_t N, typename F>
fourvec OniumRecoRate22<N, F>::calculate_fourvec(std::vector<double> parameters){
    return fourvec::unity();
}
template <size_t N, typename F>
tensor OniumRecoRate22<N, F>::calculate_tensor(std::vector<double> parameters){
    return tensor::unity();
}*/




template class OniumDissoRate22<2,double(*)(double, void*)>;

//template class OniumRecoRate22<3, double(*)(double, std::size_t, void*)>;
