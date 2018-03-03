#include "Xsection.h"
#include "matrix_elements.h"
#include "integrator.h"
#include "minimizer.h"
#include "sampler.h"
#include "predefine.h"
#include <fstream>
#include "random.h"
#include "approx_functions.h"

template<>
Xsection<2, double(*)(const double, void*)>::
	Xsection(std::string Name, std::string configfile, 
			double(*f)(const double, void*)):
StochasticBase<2>(Name+"/xsection", configfile), 
_f(f), fast_exp_(0., 15., 1000)
{
	// read configfile
	boost::property_tree::ptree config;
	std::ifstream input(configfile);
	read_xml(input, config);

	std::vector<std::string> strs;
	boost::split(strs, Name, boost::is_any_of("/"));
	std::string model_name = strs[0];
	std::string process_name = strs[1];
	auto tree = config.get_child(model_name+"."+process_name);
	_mass = tree.get<double>("mass");
	
	// Set Approximate function for X and dX_max
	StochasticBase<2>::_ZeroMoment->SetApproximateFunction(approx_X22);
	StochasticBase<2>::_FunctionMax->SetApproximateFunction(approx_dX22_max);
}

template<>
Xsection<3, double(*)(const double *, void*)>::
	Xsection(std::string Name, std::string configfile, 
			double(*f)(const double*, void*)):
StochasticBase<3>(Name+"/xsection", configfile), 
_f(f), fast_exp_(0., 15., 1000)
{
	// read configfile
	boost::property_tree::ptree config;
	std::ifstream input(configfile);
	read_xml(input, config);

	std::vector<std::string> strs;
	boost::split(strs, Name, boost::is_any_of("/"));
	std::string model_name = strs[0];
	std::string process_name = strs[1];
	auto tree = config.get_child(model_name+"."+process_name);
	_mass = tree.get<double>("mass");
	
	// Set Approximate function for X and dX_max
	StochasticBase<3>::_ZeroMoment->SetApproximateFunction(approx_X23);
	StochasticBase<3>::_FunctionMax->SetApproximateFunction(approx_dX23_max);
}

/*****************************************************************/
/*************************sample dX/dPS **************************/
/*****************************************************************/
/*------------------Implementation for 2 -> 2--------------------*/
template<>
void Xsection<2, double(*)(const double, void*)>::
		sample(std::vector<double> parameters, 
				std::vector< fourvec > & FS){
	double sqrts = parameters[0], temp = parameters[1];
	double s = std::pow(sqrts,2);
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
						StochasticBase<2>::GetFmax(parameters).s);
	double t = temp*temp*(1.-fast_exp_(-w));
	// sample phi
	double phi = Srandom::dist_phi(Srandom::gen);
	double cosphi = std::cos(phi), sinphi = std::sin(phi);
	double E = (s+_mass*_mass)/2./sqrts; // EQ = (s+M^2)/2sqrts
	double p = (s-_mass*_mass)/2./sqrts; // pQ = (s-M^2)/2sqrts
	double costheta = 1. + t/(2.*p*p); // deflection angle
	double sintheta = std::sqrt(1.-costheta*costheta);
	FS.resize(2); // HQ + light parton
	FS[0] = {E, p*sintheta*cosphi, p*sintheta*sinphi, p*costheta};
	FS[1] = {p, -p*sintheta*cosphi, -p*sintheta*sinphi, -p*costheta};
}

/*------------------Implementation for 2 -> 3--------------------*/
template<>
void Xsection<3, double(*)(const double*, void*)>::
	sample(std::vector<double> parameters,
			std::vector< fourvec > & FS){
	double sqrts = parameters[0], temp = parameters[1], 
		   delta_t = parameters[2];
	double s = sqrts*sqrts;
	auto dXdPS = [s, temp, delta_t, this](const double * PS){
		double M = this->_mass;
		double params[4] = {s, temp, M, delta_t};
		return this->_f(PS, params)/2./(s-_mass*_mass);
	};
	double Qmax = (s-_mass*_mass)/2./sqrts;
	double umax = std::log(1.+Qmax/temp);
	double xmin[4] = {0., -1.,  -1., 0.};
	double xmax[4] = {umax, 1., 1., 2.*M_PI};
	double fmax = StochasticBase<3>::GetFmax(parameters).s;
	auto res = sample_nd(dXdPS, 4, {{xmin[0], xmax[0]}, {xmin[1], xmax[1]}, 
									{xmin[2], xmax[2]}, {xmin[3], xmax[3]}},
									fmax);
}

/*****************************************************************/
/*******************find max of dX/dPS ***************************/
/*****************************************************************/
/*------------------Implementation for 2 -> 2--------------------*/
template<>
scalar Xsection<2, double(*)(const double, void*)>::
	find_max(std::vector<double> parameters){
    double s = std::pow(parameters[0],2), temp = parameters[1];
    // transform w = -log(1-t/Temp^2) 
    // since nearly everything happens when t ~ T^2
	auto minus_dXdw = [s, temp, this](const double w) {
        double T2 = temp*temp;
		double params[3] = {s, temp, this->_mass};
		double t = T2*(1.-fast_exp_(-w));
		double Jacobian = T2 - t;
		return -(this->_f(t, params)*Jacobian);
	};
	double tmin = -std::pow(s-_mass*_mass, 2)/s, tmax=0;
    double wmin = -std::log(1.-tmin/temp/temp), 
		   wmax = -std::log(1.-tmax/temp/temp+1e-9);
	double res = -minimize_1d(minus_dXdw, {wmin, wmax}, 0.0001, 100, 20);
	return scalar{res*1.5};
}
/*------------------Implementation for 2 -> 3--------------------*/
template<>
scalar Xsection<3, double(*)(const double*, void*)>::
	find_max(std::vector<double> parameters){
	double sqrts = parameters[0], temp = parameters[1], 
		   delta_t = parameters[2];
	
	double s = sqrts*sqrts;
	auto dXdPS = [s, temp, delta_t, this](const double * PS){
		double M = this->_mass;
		double params[4] = {s, temp, M, delta_t};
		return this->_f(PS, params)/2./(s-_mass*_mass);
	};
	auto nega_dXdPS = [s, temp, delta_t, this](const double * PS){
		double M = this->_mass;
		double params[4] = {s, temp, M, delta_t};
		return -this->_f(PS, params)/2./(s-_mass*_mass);
	};
	double Qmax = (s-_mass*_mass)/2./sqrts;
	double umax = std::log(1.+Qmax/temp);
	auto startloc = MC_maximize(dXdPS, 4, {{umax*0.4,umax*0.6}, {-0.2, 0.2}, 
									  {0.0, 0.8}, {-0.5, 0.5}}, 400);
	std::vector<double> step = {umax/20., 0.1, 0.05, 0.1};
	double L[4] = {0, -1, -1, 0};
	double H[4] = {umax, 1, 1, 2*M_PI};
	for(int i=0; i<4; i++){
		double dx = std::min(H[i]-startloc[i], startloc[i]-L[i])/2.;
		step[i] = std::min(dx, step[i]);
	}
    double val = -minimize_nd(nega_dXdPS, 4, startloc, step,
										2000, 1e-12);
	return scalar{val*1.1};
}

/*****************************************************************/
/*************************Integrate dX ***************************/
/*****************************************************************/
/*------------------Implementation for 2 -> 2--------------------*/
template<>
scalar Xsection<2, double(*)(const double, void*)>::
		calculate_scalar(std::vector<double> parameters){
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
/*------------------Implementation for 2 -> 3--------------------*/
template<>
scalar Xsection<3, double(*)(const double*, void*)>::
				calculate_scalar(std::vector<double> parameters){
	double sqrts = parameters[0], temp = parameters[1], 
		   delta_t = parameters[2];
	double s = sqrts*sqrts;
	auto dXdPS = [s, temp, delta_t, this](const double * PS){
		double M = this->_mass;
		double params[4] = {s, temp, M, delta_t};
		return this->_f(PS, params)/2./(s-_mass*_mass);
	};
	double Qmax = (s-_mass*_mass)/2./sqrts;
	double umax = std::log(1.+Qmax/temp);
	double xmin[4] = {0., -1., -1., 0.};
	double xmax[4] = {umax, 1., 1., 2.*M_PI};
	double error;
	double res = vegas(dXdPS, 4, xmin, xmax, error);
	return scalar{res};
}
/*****************************************************************/
/**************Integrate dX \Delta p^mu***************************/
/*****************************************************************/
/*------------------Default Implementation-----------------------*/
template<size_t N, typename F>
fourvec Xsection<N, F>::calculate_fourvec(std::vector<double> parameters){
	return fourvec{0., 0., 0., 0.};
}
/*------------------Implementation for 2 -> 2--------------------*/
template<>
fourvec Xsection<2, double(*)(const double, void*)>::
		calculate_fourvec(std::vector<double> parameters){
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

/*****************************************************************/
/**************Integrate dX \Delta p^mu*p^nu *********************/
/*****************************************************************/
/*------------------Default Implementation-----------------------*/
template<size_t N, typename F>
tensor Xsection<N, F>::calculate_tensor(std::vector<double> parameters){
	return tensor{0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,};
}
/*------------------Implementation for 2 -> 2--------------------*/
template<>
tensor Xsection<2, double(*)(const double, void*)>::
	calculate_tensor(std::vector<double> parameters){
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
// instance:
template class Xsection<2, double(*)(const double, void*)>;
template class Xsection<3, double(*)(const double*, void*)>;
