#ifndef MATRIX_ELEMENTS_H
#define MATRIX_ELEMENTS_H
#include <cstdlib>
#include <vector>
#include "lorentz.h"
//======Langevin section===============
extern double A, B;
double kperp(double E, double M, double T);
double kpara(double E, double M, double T);
void initialize_transport_coeff(double A, double B);
void postpoint_update( double dt, double M, double T, std::vector<double> v, const fourvec & pIn, fourvec & pOut);
void Ito_update( double dt, double M, double T, std::vector<double> v, const fourvec & pIn, fourvec & pOut);

//======running coupling=======================================================
double alpha_s(double Q2, double T);
double f_LPM(double x);
//======Debye mass class=======================================================
class Debye_mass{
private:
	const double TL, TH;
	const size_t NT;
	const double dT;
	const unsigned int type;
	double * mD2;
public:
	Debye_mass(const unsigned int _type);
	~Debye_mass(){delete[] mD2;};
	double get_mD2(double T);
};

extern Debye_mass * t_channel_mD2;
extern double renormalization_scale;

//=====For external initialization of debye mass==============================
void initialize_mD_and_scale(const unsigned int type, double scale, double alpha_s_fixed);

//========resonance?? Q+q-->D-->Q+q============
double dX_res_dt(const double t, void * params);

//=============Baisc function for g+q --> g+q==================================
double M2_gq2gq(const double t, void * params);
double dX_gq2gq_dt(const double t, void * params);

//=============Baisc function for g+g --> g+g==================================
double M2_gg2gg(const double t, void * params);
double dX_gg2gg_dt(const double t, void * params);

//=============Baisc function for Q+q --> Q+q==================================
double M2_Qq2Qq(const double t, void * params);
double M2_Qq2Qq_rad(const double t, void * params);
double dX_Qq2Qq_dt(const double t, void * params);

//=============Baisc function for Q+g --> Q+g==================================
double M2_Qg2Qg(const double t, void * params);
double M2_Qg2Qg_rad(const double t, void * params);
double dX_Qg2Qg_dt(const double t, void * params);

//=============Baisc function for Q+q --> Q+q+g==================================
double M2_Qq2Qqg(const double * x_, void * params_);
//=============Baisc function for Q+g --> Q+g+g==================================
double M2_Qg2Qgg(const double * x_, void * params_);

//=============Baisc function for Q+q+g --> Q+q==================================
double M2_Qqg2Qq(const double * x_,  void * params_);
//=============Baisc function for Q+g+g --> Q+g==================================
double M2_Qgg2Qg(const double * x_, void * params_);

// diffusion-like Q -> Q + g process
// July-06-2018 Yingru
double LGV_Q2Qg(const double * x_, void *params_);  // params_={E, T, M, delta_t}, x_={gluon_x, gluon_y}
double LGV_Qg2Q(const double * x_, void *params_);  // params_={E, T, M, delta_t}, x_={gluon_x, gluon_y}
#endif
