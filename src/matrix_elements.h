#ifndef MATRIX_ELEMENTS_H
#define MATRIX_ELEMENTS_H
#include <cstdlib>
#include <vector>
#include "lorentz.h"

////////// 1 <--> 2 //////////////////
// diffusion induced process
double LGV_q2qg(const double * x_, void *params_); 
double LGV_qg2q(const double * x_, void *params_); 
 
double LGV_g2gg(const double * x_, void *params_); 
double LGV_gg2g(const double * x_, void *params_);  
double LGV_g2qqbar(const double * x_, void *params_); 

////////// 2 <--> 2 //////////////////
double M2_gq2gq(const double t, void * params);
double dX_gq2gq(const double t, void * params);

double M2_gg2gg(const double t, void * params);
double dX_gg2gg(const double t, void * params);

double M2_qq2qq(const double t, void * params);
double dX_qq2qq(const double t, void * params);

double M2_qg2qg(const double t, void * params);
double dX_qg2qg(const double t, void * params);

double M2_qg2qg_stu(const double t, void * params);
double dX_qg2qg_stu(const double t, void * params);

double M2_gg2qqbar(const double t, void * params);
double dX_gg2qqbar(const double t, void * params);

double M2_qqbar2qqbar_diff(const double t, void * params);
double dX_qqbar2qqbar_diff(const double t, void * params);

////////// 2 --> 3 //////////////////
double M2_qq2qqg(const double * x_, void * params_);
double dX_qq2qqg(const double * x_, void * params_);

double M2_qg2qgg(const double * x_, void * params_);
double dX_qg2qgg(const double * x_, void * params_);

double M2_gq2gqg(const double * x_, void * params_);
double dX_gq2gqg(const double * x_, void * params_);

double M2_gg2ggg(const double * x_, void * params_);
double dX_gg2ggg(const double * x_, void * params_);

double M2_qg2qqqbar(const double * x_, void * params_);
double dX_qg2qqqbar(const double * x_, void * params_);


double M2_gq2qqqbar(const double * x_, void * params_);
double dX_gq2qqqbar(const double * x_, void * params_);

double M2_gg2qgqbar(const double * x_, void * params_);
double dX_gg2qgqbar(const double * x_, void * params_);

double M2_gg2gqqbar(const double * x_, void * params_);
double dX_gg2gqqbar(const double * x_, void * params_);

#endif
