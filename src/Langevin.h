#ifndef LANGEVIN_H
#define LANGEVIN_H

#include <vector>
#include "lorentz.h"

extern double A, B;
double qhat_pQCD(int pid, double E, double T);
double dqhat_pQCD_dp2(int pid, double E, double T);
//double kperp(double E, double M, double T);
//double kpara(double E, double M, double T);
void initialize_transport_coeff(double A, double B);

void Ito_update(int pid, double dt, double M, double T, std::vector<double> v, const fourvec & pIn, fourvec & pOut);
#endif
