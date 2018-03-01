#include "Rate.h"
#include "integrator.h"
#include "sampler.h"
#include "minimizer.h"

#include <random>

template <size_t N1, size_t N2, typename F>
Rate<N1, N2, F>::Rate(std::string Name, boost::property_tree::ptree config, F f):
StochasticBase<N1>(Name+"/rate", config),
X(std::make_shared<Xsection<N2, F>>(Name+"/xsection", config, f) ){
	auto tree = config.get_child(Name);
	_mass = tree.get<double>("mass");
	_degen = tree.get<double>("degeneracy");
	X->init();
	X->save("table.h5");
}

/*****************************************************************/
/*************************Sample dR ******************************/
/*****************************************************************/
/*------------------Implementation for 2 -> 2--------------------*/
template <>
void Rate<2, 2, double(*)(const double, void *)>::
		sample(std::vector<double> parameters, 
			std::vector< fourvec > & final_states){
	double E = parameters[0];
	double T = parameters[1];
	double v1 = std::sqrt(1. - std::pow(_mass/E,2));
	auto dR_dxdy = [E, T, v1, this](const double * x){
		double M = this->_mass;
		double E2 = T*(std::exp(x[0])-1.), costheta = x[1];
		if (costheta > 1. || costheta < -1.) return 0.;
		double s = 2.*E2*E*(1. - v1*costheta) + M*M;
		double sqrts = std::sqrt(s);
		double Xtot = this->X->GetZeroM({sqrts,T}).s;
		double Jacobian = E2 + T;
    	return 1./E*E2*std::exp(-E2/T)*(s-M*M)*2*Xtot/16./M_PI/M_PI*Jacobian;
	};
	auto res = sample_nd(dR_dxdy, 2, {{0., 3.}, {-1., 1.}}, 
						StochasticBase<2>::GetFmax(parameters).s);
	double E2 = T*(std::exp(res[0])-1.), 
		   costheta = res[1];
	double sintheta = std::sqrt(1. - costheta*costheta);
	double s = 2.*E2*E*(1. - v1*costheta) + _mass*_mass;
	double sqrts = std::sqrt(s);
	X->sample({sqrts, T}, final_states);

    // give incoming partilce a random phi angle
	double phi = std::rand()*2.*M_PI/RAND_MAX;
    // com velocity
    double vcom[3] = { E2*sintheta/(E2+E)*cos(phi), 
						E2*sintheta/(E2+E)*sin(phi),
						 (E2*costheta+v1*E)/(E2+E)	};
	/*
	FS now is in Z-oriented CoM frame
	1) FS.rotate_back
	2) FS.boost_back
	*/
	fourvec p1{E, 0, 0, v1*E};
	auto p1com = p1.boost_to(vcom[0], vcom[1], vcom[2]);
	for(auto & p: final_states){
		p = p.rotate_back(p1com);
		p = p.boost_back(vcom[0], vcom[1], vcom[2]);
	}
}
/*------------------Implementation for 2 -> 3--------------------*/
template <>
void Rate<3, 3, double(*)(const double*, void *)>::
		sample(std::vector<double> parameters, 
			std::vector< fourvec > & final_states){
	double E = parameters[0];
	double T = parameters[1];
	double v1 = std::sqrt(1. - std::pow(_mass/E,2));
	double delta_t = parameters[2];
	fourvec dxmu = {delta_t, 0., 0., delta_t*v1};
	double E2 = std::rand()*5./RAND_MAX*T;
	double costheta = std::rand()*2./RAND_MAX - 1.;
    double sintheta = std::sqrt(1. - costheta*costheta);
	double phi = std::rand()*2.*M_PI/RAND_MAX;
	double s = 2.*E2*E*(1. - v1*costheta) + _mass*_mass;
	double sqrts = std::sqrt(s);
    double vcom[3] = { E2*sintheta/(E2+E)*cos(phi), 
						E2*sintheta/(E2+E)*sin(phi),
						 (E2*costheta+v1*E)/(E2+E)	};
	double dt_com = (dxmu.boost_to(vcom[0], vcom[1], vcom[2])).t();
	X->sample({sqrts, T, dt_com}, final_states);
}

template <size_t N1, size_t N2, typename F>
scalar Rate<N1, N2, F>::find_max(std::vector<double> parameters){
	double E = parameters[0];
	double T = parameters[1];
	auto dR_dxdy = [E, T, this](const double * x){
		double M = this->_mass;
		double v1 = std::sqrt(1. - M*M/E/E);
		double E2 = T*(std::exp(x[0])-1.), costheta = x[1];
		if (E2 < 0. || costheta > 1. || costheta < -1.) return 0.;
		double s = 2.*E2*E*(1. - v1*costheta) + M*M;
		double sqrts = std::sqrt(s);
		double Xtot = this->X->GetZeroM({sqrts,T}).s;
		double Jacobian = E2 + T;
    	return -1./E*E2*std::exp(-E2/T)*(s-M*M)*2*Xtot/16./M_PI/M_PI*Jacobian;
	};
	// use f(E(x), y)*dE/dx, x = log(1+E/T), y = costheta
    // x start from 1, y start from 0
    // x step 0.3, cosphi step 0.3
    // save a slightly larger fmax
	auto val = -minimize_nd(dR_dxdy, 2, {1., 0.}, {0.3, 0.3})*1.1;
    return scalar{val};
}

/*****************************************************************/
/*************************Integrate dR ***************************/
/*****************************************************************/
/*------------------Implementation for 2 -> 2--------------------*/
template <>
scalar Rate<2, 2, double(*)(const double, void*)>::
		calculate_scalar(std::vector<double> parameters){
	double E = parameters[0];
	double T = parameters[1];
	auto code = [E, T, this](const double * x){
		double M = this->_mass;
		double v1 = std::sqrt(1. - M*M/E/E);
		double E2 = x[0], costheta = x[1];
		double s = 2.*E2*E*(1. - v1*costheta) + M*M;
		double sqrts = std::sqrt(s);
		double Xtot = this->X->GetZeroM({sqrts,T}).s;
    	std::vector<double> res{1./E*E2*std::exp(-E2/T)*(s-M*M)*2*Xtot/16./M_PI/M_PI};
		return res;
	};
	double xmin[2] = {0., -1.};
	double xmax[2] = {10.*T,1.};
	double err;
	auto val = quad_nd(code, 2, 1, xmin, xmax, err);
	return scalar{_degen*val[0]};
}
/*------------------Implementation for 2 -> 3--------------------*/
template <>
scalar Rate<3, 3, double(*)(const double*, void*)>::
		calculate_scalar(std::vector<double> parameters){
	double E = parameters[0];
	double T = parameters[1];
	double delta_t = parameters[2];
	auto code = [E, T, delta_t, this](const double * x){
		double M = this->_mass;
		double v1 = std::sqrt(1. - M*M/E/E);
		double E2 = x[0], costheta = x[1];
		double s = 2.*E2*E*(1. - v1*costheta) + M*M;
		double sqrts = std::sqrt(s);
		double Xtot = this->X->GetZeroM({sqrts,T, delta_t}).s;
    	std::vector<double> res{1./E*E2*std::exp(-E2/T)*(s-M*M)*2*Xtot/16./M_PI/M_PI};
		return res;
	};
	double xmin[2] = {0., -1.};
	double xmax[2] = {10.*T,1.};
	double err;
	auto val = quad_nd(code, 2, 1, xmin, xmax, err);
	return scalar{_degen*val[0]};
}

template <size_t N1, size_t N2, typename F>
fourvec Rate<N1, N2, F>::calculate_fourvec(std::vector<double> parameters){
	double E = parameters[0];
	double T = parameters[1];
	auto code = [E, T, this](const double * x){
		double M = this->_mass;
		double v1 = std::sqrt(1. - M*M/E/E);
		double E2 = x[0], costheta = x[1];
		double s = 2.*E2*E*(1. - v1*costheta) + M*M;
		double sintheta = std::sqrt(1. - costheta*costheta);
		double sqrts = std::sqrt(s);
		double vcom[3] = {E2*sintheta/(E2+E), 0., (E2*costheta+v1*E)/(E2+E)};
		fourvec p1{E, 0, 0, v1*E};
		// A vector in p1z(com)-oriented com frame
		auto fmu0 = this->X->GetFirstM({sqrts,T}); 
		// rotate it back from p1z(com)-oriented com frame
		auto fmu1 = fmu0.rotate_back(p1.boost_to(vcom[0], vcom[1], vcom[2]));
		// boost back to the matter frame		
		auto fmu2 = fmu1.boost_back(vcom[0], vcom[1], vcom[2]);
		double common = 1./E*E2*std::exp(-E2/T)*(s-M*M)*2./16./M_PI/M_PI;
		fmu2 = fmu2 * common;
		// Set tranverse to zero due to azimuthal symmetry;
		std::vector<double> res{fmu2.t(), fmu2.z()};
		return res;
	};
	double xmin[2] = {0., -1.};
	double xmax[2] = {5.*T, 1.};
	double err;
	auto val = quad_nd(code, 2, 2, xmin, xmax, err);
	return fourvec{_degen*val[0], 0.0, 0.0, _degen*val[1]};
}

template <size_t N1, size_t N2, typename F>
tensor Rate<N1, N2, F>::calculate_tensor(std::vector<double> parameters){
	double E = parameters[0];
	double T = parameters[1];
	auto code = [E, T, this](const double * x){
		double M = this->_mass;
		double v1 = std::sqrt(1. - M*M/E/E);
		double E2 = x[0], costheta = x[1], phi = x[2];
		double s = 2.*E2*E*(1. - v1*costheta) + M*M;
		double sintheta = std::sqrt(1. - costheta*costheta);
		double sqrts = std::sqrt(s);
		double vcom[3] = {E2*sintheta/(E2+E)*cos(phi), E2*sintheta/(E2+E)*sin(phi),
						 (E2*costheta+v1*E)/(E2+E)};
		fourvec p1{E, 0, 0, v1*E};
		// A vector in p1z(com)-oriented com frame
		auto fmunu0 = this->X->GetSecondM({sqrts,T}); 
		// rotate it back from p1z(com)-oriented com frame
		auto fmunu1 = fmunu0.rotate_back(p1.boost_to(vcom[0], vcom[1], vcom[2]));
		// boost back to the matter frame		
		auto fmunu2 = fmunu1.boost_back(vcom[0], vcom[1], vcom[2]);
		double common = 1./E*E2*std::exp(-E2/T)*(s-M*M)*2./32./std::pow(M_PI, 3);
		fmunu2 = fmunu2 * common;
		// Set tranverse to zero due to azimuthal symmetry;
		std::vector<double> res{fmunu2.T[0][0], fmunu2.T[1][1],
								fmunu2.T[2][2], fmunu2.T[3][3]};
		return res;
	};
	double xmin[3] = {0., -1., -M_PI};
	double xmax[3] = {5.*T, 1., M_PI};
	double err;
	auto val = quad_nd(code, 3, 4, xmin, xmax, err);
	return tensor{_degen*val[0], 0., 0., 0.,
				  0., _degen*val[1], 0., 0.,
				  0., 0., _degen*val[2], 0.,
				  0., 0., 0., _degen*val[3]};
}


template class Rate<2,2,double(*)(const double, void*)>; // For 2->2 
template class Rate<3,3,double(*)(const double*, void*)>; // For 2->3 
