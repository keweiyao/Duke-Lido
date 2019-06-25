#include "Onium_Rate.h"
#include "integrator.h"
#include "sampler.h"
#include "minimizer.h"
#include "random.h"
#include "approx_functions.h"
#include "matrix_elements.h"
#include "Langevin.h"
#include <iostream>


// Quarkonium 2->2 disso rate: QQbar[nl] + g --> Q + Qbar
template<>
OniumRate22<2, double(*)(double, void *)>::
    OniumRate22(std::string Name, std::string configfile, int n, int l, double(*f)(double, void *)):
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
void OniumRate22<2, double(*)(double, void *)>::
    sample(std::vector<double> parameters, std::vector< fourvec > & final_states){
}

// Find the max of dR/dq
template <>
scalar OniumRate22<2, double(*)(double, void*)>::
		find_max(std::vector<double> parameters){
	return scalar::unity();
}


// Calculate R = int_E[nl]^inf dR/dq * dq
template <>
scalar OniumRate22<2, double(*)(double, void*)>::
		calculate_scalar(std::vector<double> parameters){
	double v = parameters[0];
	double T = parameters[1];

	// dR/dx/dy to be intergated
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
fourvec OniumRate22<N, F>::calculate_fourvec(std::vector<double> parameters){
	return fourvec::unity();
}
template <size_t N, typename F>
tensor OniumRate22<N, F>::calculate_tensor(std::vector<double> parameters){
	return tensor::unity();
}



template class OniumRate22<2,double(*)(double, void*)>;
