#ifndef MATRIX_ELEMENTS_H
#define MATRIX_ELEMENTS_H
#include <cstdlib>
#include <vector>
#include "lorentz.h"

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
