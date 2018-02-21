#include "Xsection.h"
#include "matrix_elements.h"
#include "integrator.h"
#include <boost/algorithm/string.hpp>

template<size_t N, typename F>
Xsection<N, F>::Xsection(std::string Name, 
					boost::property_tree::ptree config, F f):
StochasticBase<N>(Name, config), _f(f)
{
	std::vector<std::string> strs;
	boost::split(strs, Name, boost::is_any_of("/"));
	std::string process_name = strs[0];
	auto tree = config.get_child(process_name );
	_mass = tree.get<double>("mass");
}

template<size_t N, typename F>
void Xsection<N, F>::sample(std::vector<double> arg, 
						std::vector< std::vector<double> > & FS){
	std::cout << "sampling table" << std::endl;	
}

template<size_t N, typename F>
scalar Xsection<N, F>::calculate_scalar(std::vector<double> parameters){
	double s = std::pow(parameters[0],2), temp = parameters[1];
	double * params = new double[3];
	params[0] = s; params[1] = temp; params[2]=_mass;
	auto dXdt = [params, this](double t) {return this->_f(t, params);};
	double error, tmin = -std::pow(s-_mass*_mass, 2)/s, tmax=0.;
	double res = quad_1d(dXdt, {tmin, tmax}, error);
	delete[] params;
	return scalar{res};
}

template<size_t N, typename F>
fourvec Xsection<N, F>::calculate_fourvec(std::vector<double> parameters){
	double sqrts = parameters[0], temp = parameters[1];
	double s = std::pow(sqrts,2);
	double * params = new double[3];
	params[0] = s; params[1] = temp; params[2]=_mass;
	double error, tmin = -std::pow(s-_mass*_mass, 2)/s, tmax=0.;
	const double p0 = (s-_mass*_mass)/2./sqrts;
	// <dpz*X>
	auto dpz_dXdt = [params, p0, this](double t) {
		double cos_theta13 = 1. + t/(2*p0*p0);
		double dpz = p0*(cos_theta13-1.);
		return this->_f(t, params)*dpz;
	};
	double dpz = quad_1d(dpz_dXdt, {tmin, tmax}, error);
	delete[] params;
	return fourvec{0., 0., 0., dpz};
}

template<size_t N, typename F>
tensor Xsection<N, F>::calculate_tensor(std::vector<double> parameters){
	double sqrts = parameters[0], temp = parameters[1];
	double s = std::pow(sqrts,2);
	double * params = new double[3];
	params[0] = s; params[1] = temp; params[2]=_mass;
	double error, tmin = -std::pow(s-_mass*_mass, 2)/s, tmax=0.;
	const double p0 = (s-_mass*_mass)/2./sqrts;
	// <dpz^2*X>
	auto dpzdpz_dXdt = [params, p0, this](double t) {
		double cos_theta13 = 1. + t/(2*p0*p0);
		double dpz = p0*(cos_theta13-1.);
		return this->_f(t, params)*dpz*dpz;
	};
	double dpzdpz = quad_1d(dpzdpz_dXdt, {tmin, tmax}, error);
	// <dpt^2*X>
	auto dptdpt_dXdt = [params, p0, this](double t) {
		double cos_theta13 = 1. + t/(2*p0*p0);
		double dptdpt = p0*p0*(1.-cos_theta13*cos_theta13);
		return this->_f(t, params)*dptdpt;
	};
	double dptdpt = quad_1d(dptdpt_dXdt, {tmin, tmax}, error);
	delete[] params;
	return tensor{0., 	0., 		0., 		0.,
				  0., 	dptdpt/2., 	0., 		0.,
				  0., 	0., 		dptdpt/2.,	0.,
				  0., 	0., 		0., 		dpzdpz};
}

template class Xsection<2, double(*)(double, void*)>;
