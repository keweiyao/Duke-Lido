#ifndef MATRIX_ELEMENTS_H
#define MATRIX_ELEMENTS_H
#include <cstdlib>
#include <vector>
#include "lorentz.h"


// q + (q') --> q + q'       ---> t-channel
// q + (q) --> q + q         ---> t-channel
// q + (qbar) --> q + qbar   ---> t-channel
// q + (g) --> q + g         ---> t-channel

// g + (g) --> q + qbar        ---> s/u-channel
// q + (qbar) --> g + g 		 ---> s/u-channel
// q + (qbar) --> q' + qbar'  ---> s/u-channel
// q' + (qbar') --> q + qbar  ---> s/u-channel

// q + (q) --> q + q + g
// q + (g) --> q + g + g
// g + (g) --> g + g + g
// q --> q + g

// g --> q + qbar
// q + (qbar) --> g

// g + (q or qbar) --> g + q
double M2_gq2gq(const double t, void * params);
double dX_gq2gq_dt(const double t, void * params);

// g + (g) --> g + g
double M2_gg2gg(const double t, void * params);
double dX_gg2gg_dt(const double t, void * params);

// Q + (q or qbar) --> Q + q
double M2_Qq2Qq(const double t, void * params);
double dX_Qq2Qq_dt(const double t, void * params);

// Q + (q or qbar) --> Q + q, t-channel only
double M2_Qq2Qq_rad(const double t, void * params);

// Q + (g) --> Q + g
double M2_Qg2Qg(const double t, void * params);
double dX_Qg2Qg_dt(const double t, void * params);

// Q + (g) --> Q + g, t-channel only
double M2_Qg2Qg_rad(const double t, void * params);

// Q + (q or qbar) --> Q + q +g
double M2_Qq2Qqg(const double * x_, void * params_);

// Q + (g) --> Q + g + g
double M2_Qg2Qgg(const double * x_, void * params_);

// Q + (q or qbar) + (g) --> Q + q
double M2_Qqg2Qq(const double * x_,  void * params_);

// Q + (g) + (g) --> Q + g
double M2_Qgg2Qg(const double * x_, void * params_);

// Q --> Q + g
double LGV_Q2Qg(const double * x_, void *params_);

// Q + (g) --> Q
double LGV_Qg2Q(const double * x_, void *params_); 
#endif
