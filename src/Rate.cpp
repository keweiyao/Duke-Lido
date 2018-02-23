#include "Rate.h"
#include "integrator.h"

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

template <size_t N1, size_t N2, typename F>
void Rate<N1, N2, F>::sample(std::vector<double> arg, 
						std::vector< std::vector<double> > & FS){
	std::cout << "sampling table" << std::endl;	
}

template <size_t N1, size_t N2, typename F>
scalar Rate<N1, N2, F>::find_max(std::vector<double> parameters){
    return scalar{1.0};
}

template <size_t N1, size_t N2, typename F>
scalar Rate<N1, N2, F>::calculate_scalar(std::vector<double> parameters){
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


template class Rate<2,2,double(*)(double, void*)>;
