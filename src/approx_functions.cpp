#include "approx_functions.h"

// Xsection
scalar approx_X22(std::vector<double> params){
	double sqrts = params[0];
	double T = params[1];
	return scalar{1.0/T/T};
}

scalar approx_X22QQbar(std::vector<double> params){
	double lnsqrts = params[0];
        double s = std::exp(2*lnsqrts);
	return scalar{std::log(1+s/6.76)/s};
}

scalar approx_dX22_max(std::vector<double> params){
	double T = params[1];
	return scalar{1.0/std::pow(T, 2)};
}
