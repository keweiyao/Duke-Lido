#ifndef APPROX_FUNCTIONS_H
#define APPROX_FUNCTIONS_H
#include <vector>
#include "lorentz.h"
scalar approx_X22(std::vector<double> params){
	double sqrts = params[0];
	double T = params[1];
	return scalar{1.0/T/T};
}

scalar approx_dX22_max(std::vector<double> params){
	double sqrts = params[0];
	double T = params[1];
	return scalar{1.0/std::pow(T, 2)};
}

scalar approx_X23(std::vector<double> params){
	double sqrts = params[0];
	double T = params[1];
	double delta_t = params[2];
	return scalar{1./(1.+1./std::pow(3*delta_t*T,2))
					/std::pow(T, 2)*sqrts*sqrts};
}

scalar approx_dX23_max(std::vector<double> params){
	double sqrts = params[0];
	double T = params[1];
	double delta_t = params[2];
	return scalar{1./(1.+1./std::pow(delta_t*T,2))
					/std::pow(T, 3)*std::pow(sqrts*sqrts, 2)};
}

#endif
