#ifndef APPROX_FUNCTIONS_H
#define APPROX_FUNCTIONS_H
#include<vector>

double approx_X22(std::vector<double> params){
	double sqrts = params[0];
	double T = params[1];
	double M = params[2];
	return 1.0/(sqrts*sqrts-M*M)/T/T;
}

double approx_dX22_max(std::vector<double> params){
	double sqrts = params[0];
	double T = params[1];
	double M = params[2];
	return 1.0/(sqrts*sqrts-M*M)/std::pow(T, 4);
}

double approx_dX23_max(std::vector<double> params){
	double sqrts = params[0];
	double T = params[1];
	double M = params[2];
	return 1./(1.+1./std::pow(delta_t*T,2))
					/std::pow(T, 2)*std::pow(1.-std::pow(M/sqrts,2), 2);
}

double approx_X23(std::vector<double> params){
	double sqrts = params[0];
	double T = params[1];
	double M = params[2];
	double delta_t = params[3]
	return 1./(1.+1./std::pow(3*delta_t*T,2))
					/std::pow(T, 6)*std::pow(1.-std::pow(M/sqrts,2), 2);
}

#endif
