#include <cmath>
#include <iostream>
#include "predefine.h"
#include "matrix_elements.h"
#include "lorentz.h"
#include "simpleLogger.h"

///	g+g --> g+g
double M2_gg2gg(const double t, void * params){
	// unpacking parameters
	double * p = static_cast<double*>(params);
	double s = p[0], Temp = p[1], M2 = 0.;
	// Deybe mass (t-channel)
	double mt2 = t_channel_mD2->get_mD2(Temp);
	// define energy scales for each channel
	double Q2s = s - M2, Q2t = std::min(t, -cut*mt2), Q2u = M2 - s - t;
	// define coupling constant for each channel
	double At = alpha_s(Q2t, Temp);

	double result = 72.*M_PI*M_PI*At*At*(-Q2s*Q2u/Q2t/Q2t);

	if (result < 0.) return 0.;
	return result;
}
double dX_gg2gg_dt(const double t, void * params){
	double * p = static_cast<double*>(params);
	double s = p[0], M2 = 0.;
	return M2_gg2gg(t, params)/c16pi/(std::pow(s-M2, 2));
}
///	g+q --> g+q
double M2_gq2gq(const double t, void * params){
	// unpacking parameters
	double * p = static_cast<double*>(params);
	double s = p[0], Temp = p[1], M2 = 0.;
	// Deybe mass (t-channel)
	double mt2 = t_channel_mD2->get_mD2(Temp);
	// define energy scales for each channel
	double Q2s = s - M2, Q2t = std::min(t, -cut*mt2), Q2u = M2 - s - t;
	// define coupling constant for each channel
	double At = alpha_s(Q2t, Temp);

	double result = At*At*64.*M_PI*M_PI/9.*(Q2s*Q2s+Q2u*Q2u)*( 9./4./Q2t/Q2t);
	if (result < 0.) return 0.;
	return result;
}
double dX_gq2gq_dt(const double t, void * params){
	double * p = static_cast<double*>(params);
	double s = p[0], M2 = 0.;
	return M2_gq2gq(t, params)/c16pi/(std::pow(s-M2, 2));
}


/// Q+q --> Q+q
double M2_Qq2Qq(const double t, void * params){
	// unpacking parameters
	double * p = static_cast<double*>(params);
	double s = p[0], Temp = p[1], M2 = p[2]*p[2];
	// Deybe mass (t-channel)
	double mt2 = t_channel_mD2->get_mD2(Temp);
	// define energy scales for each channel
	double Q2s = s - M2, Q2t = std::min(t, -cut*mt2), Q2u = M2 - s - t;
	double At = alpha_s(Q2t, Temp);

	double result = c64d9pi2*At*At*(Q2u*Q2u + Q2s*Q2s + 2.*M2*t)/Q2t/Q2t;
	if (result < 0.) return 0.;
	else return result;
}

/// Q+q --> Q+q for radiation process, rightnow, this is the same as Q+q-->Q+q,
/// exept that we don't restrict (-t) > mD2 here, instead, we require later
/// that the emitted gluon has an energy larger than mD, excludes (-t) < epsilon
double M2_Qq2Qq_rad(const double t, void * params){
	// unpacking parameters
	double * p = static_cast<double*>(params);
	double s = p[0], Temp = p[1], M2 = p[2]*p[2];
	// Deybe mass (t-channel)
	double mt2 = t_channel_mD2->get_mD2(Temp);
	// define energy scales for each channel
	double Q2s = s - M2, Q2t = std::min(t, -cut*mt2), Q2u = M2 - s - t;
	double At = alpha_s(Q2t, Temp);

	double result = c64d9pi2*At*At*(Q2u*Q2u + Q2s*Q2s + 2.*M2*t)/Q2t/Q2t;

	if (result < 0.) return 0.;
	else return result;
}

double dX_Qq2Qq_dt(const double t, void * params){
	double * p = static_cast<double*>(params);
	double s = p[0], M2 = p[2]*p[2];
	return M2_Qq2Qq(t, params)/c16pi/std::pow(s-M2, 2);
}

///	Q+g --> Q+g
double M2_Qg2Qg(const double t, void * params){
	// unpacking parameters
	double * p = static_cast<double*>(params);
	double s = p[0], Temp = p[1], M2 = p[2]*p[2];
	// Deybe mass (t-channel)
	double mt2 = t_channel_mD2->get_mD2(Temp);
	// define energy scales for each channel
	double Q2s = s - M2, Q2t = std::min(t, -cut*mt2), Q2u = M2 - s - t;
	// define coupling constant for each channel
	double At = alpha_s(Q2t, Temp),
		   Au = alpha_s(Q2u, Temp),
		   As = alpha_s(Q2s, Temp);
	double Q2s_reg = Q2s + mt2;
	double Q2u_reg = Q2u>0?(Q2u + mt2):(Q2u-mt2);

	double result = 0.0;
	// t*t
	result += 2.*At*At * Q2s*(-Q2u)/Q2t/Q2t;
	// s*s
	result += c4d9*As*As *
			( Q2s*(-Q2u) + 2.*M2*(Q2s + 2.*M2) ) / std::pow(Q2s_reg, 2);
	// u*u
	result += c4d9*Au*Au *
			( Q2s*(-Q2u) + 2.*M2*(Q2u + 2.*M2) ) / std::pow(Q2u_reg, 2);
	// s*u
	result += c1d9*As*Au * M2*(4.*M2 - t) / Q2s_reg / (-Q2u_reg);
	// t*s
	result += At*As * ( Q2s*(-Q2u) + M2*(Q2s - Q2u) ) / Q2t / Q2s_reg;
    // t*u
	result += -At*Au * ( Q2s*(-Q2u) - M2*(Q2s - Q2u) ) / Q2t / (-Q2u_reg);
	if (result < 0.) return 0.;
	return result*c16pi2;
}


/// Q+g --> Q+g for radiation process, rightnow, this only inclues t-channel
/// We don't restrict (-t) > mD2 here, instead, we require later that the
/// emitted gluon has an energy larger than mD, which exclude (-t) < epsilon
double M2_Qg2Qg_rad(const double t, void * params) {
	// unpacking parameters
	double * p = static_cast<double *>(params);
	double s = p[0], Temp = p[1], M2 = p[2]*p[2];
	// Deybe mass (t-channel)
	double mt2 = t_channel_mD2->get_mD2(Temp);
	double Q2s = s - M2, Q2t = std::min(t, -cut*mt2), Q2u = M2 - s - t;
	double At = alpha_s(Q2t, Temp);

	double result = c16pi2*2.*At*At * Q2s*(-Q2u)/Q2t/Q2t;

	if (result < 0.) return 0.;
	else return result;
}

double dX_Qg2Qg_dt(const double t, void * params){
	double * p = static_cast<double*>(params);
	double s = p[0], M2 = p[2]*p[2];
	return M2_Qg2Qg(t, params)/c16pi/std::pow(s-M2, 2);
}

/// Q + q --> Q + q + g
double M2_Qq2Qqg(const double * x_, void * params_){
	// unpack variables, parameters and check integration range
	double * params = static_cast<double*>(params_);
	double s = params[0];
	double sqrts = std::sqrt(s);
	double T = params[1];
	double M2 = params[2]*params[2];
	double Qmax = (s-M2)/2./sqrts;
	double log_1_ktT = x_[0], y_norm = x_[1], // in 3-body com frame---(1)
		   costheta34 = x_[2], phi34 =x_[3]; // in p3-p4 com frame---(2)
		   // (1) and (2) share the same z-direction, and are related by a boost
	// check bounds
	double kt = T*(std::exp(log_1_ktT)-1.);
	if (std::abs(costheta34)>=1.||phi34<=0.||phi34>=2.*M_PI||kt>=Qmax||kt<=0.)
		return 0.;
	double ymax = std::acosh(Qmax/kt);
	double y = y_norm*ymax;
	// construct k^mu
	fourvec kmu{kt*std::cosh(y), kt, 0, kt*std::sinh(y)};
	double s34 = s - 2.*sqrts*kmu.t();
	double sqrts34 = std::sqrt(s34);
	double Q34 = (s34-M2)/2./sqrts34, E34 = (s34+M2)/2./sqrts34;
	// construct p3, p4 in (2) frame
	double cosphi34 = std::cos(phi34), sinphi34 = std::sin(phi34),
		   sintheta34 = std::sqrt(1.-std::pow(costheta34,2));
	fourvec p3mu{E34, Q34*sintheta34*cosphi34,
				Q34*sintheta34*sinphi34, Q34*costheta34};
	fourvec p4mu{Q34, -Q34*sintheta34*cosphi34,
				-Q34*sintheta34*sinphi34, -Q34*costheta34};
	// boost p3, p4 back to 3-body CoM Frame
	double V0 = sqrts - kmu.t();
	double v34[3] = {-kmu.x()/V0, -kmu.y()/V0, -kmu.z()/V0};
	p3mu = p3mu.boost_back(v34[0], v34[1], v34[2]);
	p4mu = p4mu.boost_back(v34[0], v34[1], v34[2]);

	// q-perp-vec, q = p2-p4, qperp = -p4perp
	double qx = -p4mu.x(), qy = -p4mu.y();
	double qt2 = qx*qx + qy*qy;

	double mD2 = t_channel_mD2->get_mD2(T);
	double kt2 = kt*kt;
	double kt_qt2 = kt2 - 2.*qx*kmu.x() + qt2;
	double x = (kmu.t()+kmu.z())/sqrts,
	       xbar = (kmu.t()+std::abs(kmu.z()))/sqrts;
	double one_minus_xbar = 1.-xbar;
	
	double iD1 = 1./(kt2 + (1.-xbar)*mD2/2.),
	       iD2 = 1./(kt_qt2 + (1.-xbar)*mD2/2.);

	double t = -2.*Qmax*(p4mu.t()+p4mu.z());
	double M2_elastic = M2_Qq2Qq_rad(t, params);

	double alpha_rad = alpha_s(kt2, T);
	double Pg = alpha_rad*std::pow(one_minus_xbar, 2)
      *(kt2*std::pow(iD1-iD2, 2.) + qt2*std::pow(iD2,2) + 2.*kmu.x()*qx*iD2*(iD1-iD2));
	// Jacobian
	double J = (1.0 - M2/s34) * M_PI/8./std::pow(2*M_PI,5) * kt * (kt + T) * ymax;
	// 2->3 = 2->2 * 1->2
	return c48pi*M2_elastic*Pg*J;
}

/// Q + g --> Q + g + g
double M2_Qg2Qgg(const double * x_, void * params_){
	// unpack variables, parameters and check integration range
	double * params = static_cast<double*>(params_);
	double s = params[0];
	double sqrts = std::sqrt(s);
	double T = params[1];
	double M2 = params[2]*params[2];
	double Qmax = (s-M2)/2./sqrts;
	double log_1_ktT = x_[0], y_norm = x_[1], // in 3-body com frame---(1)
		   costheta34 = x_[2], phi34 =x_[3]; // in p3-p4 com frame---(2)
		   // (1) and (2) share the same z-direction, and are related by a boost
	// check bounds
	double kt = T*(std::exp(log_1_ktT)-1.);
	if (std::abs(costheta34)>=1.||phi34<=0.||phi34>=2.*M_PI||kt>=Qmax||kt<=0.)
		return 0.;
	double ymax = std::acosh(Qmax/kt);
	double y = y_norm*ymax;
	// construct k^mu
	fourvec kmu{kt*std::cosh(y), kt, 0, kt*std::sinh(y)};
	double s34 = s - 2.*sqrts*kmu.t();
	double sqrts34 = std::sqrt(s34);
	double Q34 = (s34-M2)/2./sqrts34, E34 = (s34+M2)/2./sqrts34;
	// construct p3, p4 in (2) frame
	double cosphi34 = std::cos(phi34), sinphi34 = std::sin(phi34),
		   sintheta34 = std::sqrt(1.-std::pow(costheta34,2));
	fourvec p3mu{E34, Q34*sintheta34*cosphi34,
				Q34*sintheta34*sinphi34, Q34*costheta34};
	fourvec p4mu{Q34, -Q34*sintheta34*cosphi34,
				-Q34*sintheta34*sinphi34, -Q34*costheta34};
	// boost p3, p4 back to 3-body CoM Frame
	double V0 = sqrts - kmu.t();
	double v34[3] = {-kmu.x()/V0, -kmu.y()/V0, -kmu.z()/V0};
	p3mu = p3mu.boost_back(v34[0], v34[1], v34[2]);
	p4mu = p4mu.boost_back(v34[0], v34[1], v34[2]);

	// q-perp-vec, q = p2-p4, qperp = -p4perp
	double qx = -p4mu.x(), qy = -p4mu.y();
	double qt2 = qx*qx + qy*qy;

	double mD2 = t_channel_mD2->get_mD2(T);
    double kt2 = kt*kt;
    double kt_qt2 = kt2 - 2.*qx*kmu.x() + qt2;
    double x = (kmu.t()+kmu.z())/sqrts,
           xbar = (kmu.t()+std::abs(kmu.z()))/sqrts;
    double one_minus_xbar = 1.-xbar;
    double iD1 = 1./(kt2 + (1.-xbar)*mD2/2.),
           iD2 = 1./(kt_qt2  + (1.-xbar)*mD2/2.);

    double t = -2.*Qmax*(p4mu.t()+p4mu.z());
    double M2_elastic = M2_Qg2Qg_rad(t, params);

    double alpha_rad = alpha_s(kt2, T);
    double Pg = alpha_rad*std::pow(one_minus_xbar, 2)
  *(kt2*std::pow(iD1-iD2, 2.) + qt2*std::pow(iD2,2) + 2.*kmu.x()*qx*iD2*(iD1-iD2));

	// Jacobian
	double J = (1.0 - M2/s34) * M_PI/8./std::pow(2*M_PI,5) * kt * (kt + T) * ymax;
	// 2->3 = 2->2 * 1->2
	return c48pi*M2_elastic*Pg*J;
}


//=============Basic for 3->2===========================================
// Not CoM frame of  1,2,k, but CoM frame of 1 and 2
// sampled final states within 3+4 CoM frame
double M2_Qqg2Qq(const double * x_, void * params_){
	// unpack variables costheta42 = x_[0]
	double costheta34 = x_[0], phi34 = x_[1];

	if (costheta34<=-1. || costheta34>=1. || phi34 <=0. || phi34 >=2.*M_PI) return 0.;
	double sintheta34 = std::sqrt(1. - costheta34*costheta34),
		sinphi34 = std::sin(phi34), cosphi34 = std::cos(phi34);
	// unpack parameters
	double * params = static_cast<double*>(params_); // s12, T, M, k, costhetak
	double s = params[0];
	double T = params[1];
	double M = params[2];
	double M2 = M*M;
	double s12 = params[3]*(s-M2) + M2; // 0 < params[3] = xinel = (s12-M2)/(s-M2) < 1
	double s1k = (params[4]*(1-s12/s)*(1-M2/s12)+M2/s12)*s; // 0 < params[4] = yinel = (s1k/s-M2/s12)/(1-s12/s)/(1-M2/s12) < 1
	double sqrts12 = std::sqrt(s12);
	double sqrts = std::sqrt(s);
	double E1 = (s12+M2)/2./sqrts12,  p1 = (s12-M2)/2./sqrts12;
	double E3 = (s+M2)/2./sqrts, p3 = (s-M2)/2./sqrts;
	double k = (s-s12)/2./sqrts12;
	double costhetak = (M2 + 2.*E1*k - s1k)/2./p1/k;
	double sinthetak = std::sqrt(1. - costhetak*costhetak);
	double kt = k*sinthetak;
	double kt2 = kt*kt;
	double kz = k*costhetak;
	double x = (k+kz)/(sqrts12+k+kz), xbar = (k+std::abs(kz))/(sqrts12+k+std::abs(kz));
  double mD2 = t_channel_mD2->get_mD2(T);

	// get final state
	fourvec Ptot{sqrts12+k, kt, 0., kz};
	double vcom[3] = {Ptot.x()/Ptot.t(), Ptot.y()/Ptot.t(),Ptot.z()/Ptot.t()};
	// final state in 34-com frame
	fourvec p3mu{E3, p3*sintheta34*cosphi34, p3*sintheta34*sinphi34, p3*costheta34};
	fourvec p4mu{p3, -p3*sintheta34*cosphi34, -p3*sintheta34*sinphi34, -p3*costheta34};

	// boost final state back to 12-com frame
	p3mu = p3mu.boost_back(vcom[0], vcom[1], vcom[2]);
	p4mu = p4mu.boost_back(vcom[0], vcom[1], vcom[2]);

	// 2->2 part
	fourvec qmu{p1-p4mu.t(), -p4mu.x(), -p4mu.y(), -p1-p4mu.z()};
	double t = -2.*p1*(p4mu.t()+p4mu.z());

	double qt2 = std::pow(qmu.x(),2) + std::pow(qmu.y(),2);

	double new_params[3] = {s, T, M};
	double M2_elastic = M2_Qq2Qq_rad(t, params);
	double x2M2 = x*x*M2;

	// 1->2
	double iD1 = 1./(kt2 + (1.-xbar)*mD2/2.);
	double iD2 = 1./(kt2 + qt2 + 2.*kt*qmu.x() + (1.-xbar)*mD2/2.);
	double Pg = 48.*M_PI*alpha_s(kt2, T)*std::pow(1.-xbar, 2)*(
			kt2*std::pow(iD1-iD2, 2) + qt2*std::pow(iD2, 2) - 2.*kt*qmu.x()*(iD1-iD2)*iD2
		);

	double detail_balance_factor = 1./16.;
	double Jacobian = 1./std::pow(2*M_PI, 2)/8.*(1. - M2/s);
	// 2->3 = 2->2 * 1->2
	return M2_elastic * Pg * Jacobian * detail_balance_factor;
}

double M2_Qgg2Qg(const double * x_, void * params_){
	// unpack variables costheta42 = x_[0]
	double costheta34 = x_[0], phi34 = x_[1];

	if (costheta34<=-1. || costheta34>=1. || phi34 <=0. || phi34 >=2.*M_PI) return 0.;
	double sintheta34 = std::sqrt(1. - costheta34*costheta34),
		sinphi34 = std::sin(phi34), cosphi34 = std::cos(phi34);
	// unpack parameters
	double * params = static_cast<double*>(params_); // s12, T, M, k, costhetak
	double s = params[0];
	double T = params[1];
	double M = params[2];
	double M2 = M*M;
	double s12 = params[3]*(s-M2) + M2; // 0 < params[3] = xinel = (s12-M2)/(s-M2) < 1
	double s1k = (params[4]*(1-s12/s)*(1-M2/s12)+M2/s12)*s; // 0 < params[4] = yinel = (s1k/s-M2/s12)/(1-s12/s)/(1-M2/s12) < 1
	double sqrts12 = std::sqrt(s12);
	double sqrts = std::sqrt(s);
	double E1 = (s12+M2)/2./sqrts12,  p1 = (s12-M2)/2./sqrts12;
	double E3 = (s+M2)/2./sqrts, p3 = (s-M2)/2./sqrts;
	double k = (s-s12)/2./sqrts12;
	double costhetak = (M2 + 2.*E1*k - s1k)/2./p1/k;
	double sinthetak = std::sqrt(1. - costhetak*costhetak);
	double kt = k*sinthetak;
	double kt2 = kt*kt;
	double kz = k*costhetak;
	double x = (k+kz)/(sqrts12+k+kz), xbar = (k+std::abs(kz))/(sqrts12+k+std::abs(kz));
    double mD2 = t_channel_mD2->get_mD2(T);

	// get final state
	fourvec Ptot{sqrts12+k, kt, 0., kz};
	double vcom[3] = {Ptot.x()/Ptot.t(), Ptot.y()/Ptot.t(),Ptot.z()/Ptot.t()};
	// final state in 34-com frame
	fourvec p3mu{E3, p3*sintheta34*cosphi34, p3*sintheta34*sinphi34, p3*costheta34};
	fourvec p4mu{p3, -p3*sintheta34*cosphi34, -p3*sintheta34*sinphi34, -p3*costheta34};

	// boost final state back to 12-com frame
	p3mu = p3mu.boost_back(vcom[0], vcom[1], vcom[2]);
	p4mu = p4mu.boost_back(vcom[0], vcom[1], vcom[2]);

	// 2->2 part
	fourvec qmu{p1-p4mu.t(), -p4mu.x(), -p4mu.y(), -p1-p4mu.z()};
	double t = -2.*p1*(p4mu.t()+p4mu.z());

	double qt2 = std::pow(qmu.x(),2) + std::pow(qmu.y(),2);

	double new_params[3] = {s, T, M};
	double M2_elastic = M2_Qg2Qg_rad(t, params);
	double x2M2 = x*x*M2;

	// 1->2
	double iD1 = 1./(kt2 + (1.-xbar)*mD2/2.);
	double iD2 = 1./(kt2 + qt2 + 2.*kt*qmu.x() + (1.-xbar)*mD2/2.);
	double Pg = 48.*M_PI*alpha_s(kt2, T)*std::pow(1.-xbar, 2)*(
			kt2*std::pow(iD1-iD2, 2) + qt2*std::pow(iD2, 2) - 2.*kt*qmu.x()*(iD1-iD2)*iD2
		);

	double detail_balance_factor = 1./16.;
	double Jacobian = 1./std::pow(2*M_PI, 2)/8.*(1. - M2/s);
	// 2->3 = 2->2 * 1->2
	return M2_elastic * Pg * Jacobian * detail_balance_factor;
}

// July-06-2019
// Diffusion-induced Radiative process Q -> Q + g
double P_q2qg(double x){
	return CF*(1 + std::pow(1-x, 2))/x;
}
double x_P_q2qg(double x){
	return CF*(1 + std::pow(1-x, 2));
}

double LGV_Q2Qg(const double * x_, void * params_){
    double *params = static_cast<double*>(params_);
    double E = params[0];
    double T = params[1];
    double M = params[2];
	double pabs = std::sqrt(E*E-M*M);
    
	double x = x_[0]; // k/p
	double y = x_[1]; // kT/k0

    double k0 = x*pabs;
    double kT = y*k0;
	double kT2 = kT*kT;
    
	double mg2 = t_channel_mD2->get_mD2(T)/2.;
	double x0 = k0/E;
	// No dead cone
	double Jacobian = 2*k0*kT;
    double dR_dxdy = alpha_s(kT2, T)/(2.*M_PI) * P_q2qg(x0)
                     * 1./std::pow(kT2+mg2, 2)
                     * Jacobian;
    return dR_dxdy;
}


// Diffusion-induced gluon absorption process Q + g -> Q
double LGV_Qg2Q(const double * x_, void * params_){
    double *params = static_cast<double*>(params_);
    double E = params[0];
    double T = params[1];
    double M = params[2];
	double pabs = std::sqrt(E*E-M*M);
    
	double xp = x_[0];// tanh(x)
	double x = std::atanh(xp); // kz/pabs
	double y = x_[1];// kT/k

	double kz = x*pabs;
	double k0 = std::abs(kz)/std::sqrt(1.001-y*y);
	double kT = y*k0;

	double kT2 = kT*kT;
	double x0 = k0/(k0+E); // k0/(E0)

	double mg2 = t_channel_mD2->get_mD2(T)/2.;

	// No dead cone
	double Jacobian = 2*(k0+E)*kT/(1-xp*xp);
    double dR_dxdy = alpha_s(kT2, T)/(2.*M_PI) * x_P_q2qg(x0)
                     * 1./std::pow(kT2+mg2, 2)
                     * Jacobian * std::exp(-k0/T);
    return dR_dxdy;
}
