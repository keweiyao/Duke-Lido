#include "Xsection.h"
#include "matrix_elements.h"
#include "integrator.h"
#include "minimizer.h"
#include "sampler.h"
#include <boost/algorithm/string.hpp>

template<size_t N, typename F>
Xsection<N, F>::Xsection(std::string Name, 
					boost::property_tree::ptree config, F f):
StochasticBase<N>(Name, config), _f(f), fast_exp_(0., 15., 1000)
{
	std::vector<std::string> strs;
	boost::split(strs, Name, boost::is_any_of("/"));
	std::string process_name = strs[0];
	auto tree = config.get_child(process_name );
	_mass = tree.get<double>("mass");
}

template<size_t N, typename F>
void Xsection<N, F>::sample(std::vector<double> parameters, 
						std::vector< std::vector<double> > & FS){
	double s = std::pow(parameters[0],2), temp = parameters[1];
    // transform w = -log(1-t/Temp^2) 
    // since nearly everything happens when t ~ T^2
	auto dXdw = [s, temp, this](double w) {
        double T2 = temp*temp;
		double params[3] = {s, temp, this->_mass};
		double t = T2*(1.-fast_exp_(-w));
		double Jacobian = T2 - t;
		return this->_f(t, params)*Jacobian;
	};
	double tmin = -std::pow(s-_mass*_mass, 2)/s, tmax=0;
    double wmin = -std::log(1.-tmin/temp/temp), 
		   wmax = -std::log(1.-tmax/temp/temp+1e-9);
	double w = sample_1d(dXdw, {wmin, wmax}, 
						StochasticBase<N>::GetFmax(parameters).s);
	double t = temp*temp*(1.-fast_exp_(-w));
	//std::cout << t << std::endl;
}

template<size_t N, typename F>
scalar Xsection<N, F>::find_max(std::vector<double> parameters){
    double s = std::pow(parameters[0],2), temp = parameters[1];
    // transform w = -log(1-t/Temp^2) 
    // since nearly everything happens when t ~ T^2
	auto minus_dXdw = [s, temp, this](double w) {
        double T2 = temp*temp;
		double params[3] = {s, temp, this->_mass};
		double t = T2*(1.-fast_exp_(-w));
		double Jacobian = T2 - t;
		return -(this->_f(t, params)*Jacobian);
	};
	double tmin = -std::pow(s-_mass*_mass, 2)/s, tmax=0;
    double wmin = -std::log(1.-tmin/temp/temp), 
		   wmax = -std::log(1.-tmax/temp/temp+1e-9);
	double res = -minimize_1d(minus_dXdw, {wmin, wmax});
 
	return scalar{res};
}

template<size_t N, typename F>
scalar Xsection<N, F>::calculate_scalar(std::vector<double> parameters){
	double s = std::pow(parameters[0],2), temp = parameters[1];
    // transform w = -log(1-t/Temp^2) 
    // since nearly everything happens when t ~ T^2
	auto dXdw = [s, temp, this](double w) {
        double T2 = temp*temp;
		double params[3] = {s, temp, this->_mass};
		double t = T2*(1.-fast_exp_(-w));
		double Jacobian = T2 - t;
		return this->_f(t, params)*Jacobian;
	};
	double error, tmin = -std::pow(s-_mass*_mass, 2)/s, tmax=0;
    double wmin = -std::log(1.-tmin/temp/temp), 
		   wmax = -std::log(1.-tmax/temp/temp+1e-9);
	double res = quad_1d(dXdw, {wmin, wmax}, error);
	return scalar{res};
}

template<size_t N, typename F>
fourvec Xsection<N, F>::calculate_fourvec(std::vector<double> parameters){
	double sqrts = parameters[0], temp = parameters[1];
	double s = std::pow(sqrts,2);
	double p0 = (s-_mass*_mass)/2./sqrts;
	// <dpz*X>
	auto dpz_dXdw = [s, temp, p0, this](double w) {
        double T2 = temp*temp;
		double params[3] = {s, temp, this->_mass};
		double t = T2*(1.-fast_exp_(-w));
		double Jacobian = T2 - t;
        double cos_theta13 = 1. + t/(2*p0*p0);
        double dpz = p0*(cos_theta13-1.);
		return this->_f(t, params)*dpz*Jacobian;
	};
	double error, tmin = -std::pow(s-_mass*_mass, 2)/s, tmax=0;
    double wmin = -std::log(1.-tmin/temp/temp),
		   wmax = -std::log(1.-tmax/temp/temp+1e-9);
	double dpz = quad_1d(dpz_dXdw, {wmin, wmax}, error);
	return fourvec{0., 0., 0., dpz};
}

template<size_t N, typename F>
tensor Xsection<N, F>::calculate_tensor(std::vector<double> parameters){
	double sqrts = parameters[0], temp = parameters[1];
	double s = std::pow(sqrts,2);
	const double p0 = (s-_mass*_mass)/2./sqrts;
	// <dpz^2*X>
	auto dpzdpz_dXdw = [s, temp, p0, this](double w) {
        double T2 = temp*temp;
		double params[3] = {s, temp, this->_mass};
		double t = T2*(1.-fast_exp_(-w));
		double Jacobian = T2 - t;
        double cos_theta13 = 1. + t/(2*p0*p0);
        double dpz = p0*(cos_theta13-1.);
		return this->_f(t, params)*dpz*dpz*Jacobian;
	};
	// <dpt^2*X>
	auto dptdpt_dXdw = [s, temp, p0, this](double w) {
		double T2 = temp*temp;
		double params[3] = {s, temp, this->_mass};
		double t = T2*(1.-std::exp(-w));
		double Jacobian = T2 - t;
        double cos_theta13 = 1. + t/(2*p0*p0);
		double dptdpt = p0*p0*(1.-cos_theta13*cos_theta13);
		return this->_f(t, params)*dptdpt*Jacobian;
	};

	double error, tmin = -std::pow(s-_mass*_mass, 2)/s, tmax=0;
    double wmin = -std::log(1.-tmin/temp/temp), 
		   wmax = -std::log(1.-tmax/temp/temp+1e-9);
    double dpzdpz = quad_1d(dpzdpz_dXdw, {wmin, wmax}, error);
	double dptdpt = quad_1d(dptdpt_dXdw, {wmin, wmax}, error);
	return tensor{0., 	0., 		0., 		0.,
				  0., 	dptdpt/2., 	0., 		0.,
				  0., 	0., 		dptdpt/2.,	0.,
				  0., 	0., 		0., 		dpzdpz};
}

template class Xsection<2, double(*)(double, void*)>;
