#include <cmath>
#include <iostream>
#include "predefine.h"
#include "matrix_elements.h"
#include "lorentz.h"
#include "simpleLogger.h"

/////////////////////// SPLITTING ////////////////////////

// Diffusion-induced gluon radiation process q -> q + g
double LGV_q2qg(const double * x_, void * params_){
  double *params = static_cast<double*>(params_);
  double EA = params[0];
  double T = params[1];
  double mA = params[2], m1 = params[3], m2 = params[4];
  double pA = std::sqrt(EA*EA-mA*mA);
  double x = x_[0]; // k+/p+
  double sintheta = x_[1]; // sintheta = kT/kabs > 0
  double tantheta = sintheta/std::sqrt(1.-sintheta*sintheta);
  double kplus = x*(EA+pA);
  double kT = std::sqrt((kplus/sintheta+m1)*(kplus/sintheta-m1))
            - kplus/tantheta;
  double kT2 = kT*kT;
  double mg2 = (CF/CA*x*x + 1. -x)*t_channel_mD2->get_mD2(T)/2.;
  double Jacobian = 2*kT2/sintheta;
  double dR_dxdy = alpha_s(kT2, T)/(2.*M_PI)
                     * P_q2qg(x)
                     / std::pow(kT2+mg2, 2)
                     * Jacobian * (x*x*CF/CA+(1-x));
  return dR_dxdy;
}

// Diffusion-induced gluon radiation process g -> g + g
double LGV_g2gg(const double * x_, void * params_){
  double *params = static_cast<double*>(params_);
  double EA = params[0];
  double T = params[1];
  double mA = params[2], m1 = params[3], m2 = params[4];
  double pA = std::sqrt(EA*EA-mA*mA);

  double x = x_[0]; // k+/p+
  double sintheta = x_[1]; // sintheta = kT/kabs > 0
  double tantheta = sintheta/std::sqrt(1.-sintheta*sintheta);
  double kplus = x*(EA+pA);
  double kT = std::sqrt((kplus/sintheta+m1)*(kplus/sintheta-m1))
            - kplus/tantheta;
  double kT2 = kT*kT;
  double mg2 = (1.-x+x*x)*t_channel_mD2->get_mD2(T)/2.;
  if (x > 0.5) return 0.0; 
  double Jacobian = 2.*kT2/sintheta;
  double dR_dxdy = alpha_s(kT2, T)/(2.*M_PI) 
                     * P_g2gg(x)
                     / std::pow(kT2+mg2, 2)
                     * Jacobian * (1.-x+x*x);
  return dR_dxdy;
}

// Diffusion-induced gluon splitting process g -> q + qbar
double LGV_g2qqbar(const double * x_, void * params_){
  double *params = static_cast<double*>(params_);
  double EA = params[0];
  double T = params[1];
  double mA = params[2], m1 = params[3], m2 = params[4];
  double pA = std::sqrt(EA*EA-mA*mA);

  double x = x_[0]; // k+/p+
  double sintheta = x_[1]; // sintheta = kT/kabs > 0
  double tantheta = sintheta/std::sqrt(1.-sintheta*sintheta);
  double kplus = x*(EA+pA);
  double kT = std::sqrt((kplus/sintheta+m1)*(kplus/sintheta-m1))
            - kplus/tantheta;
  double kT2 = kT*kT;
  double mg2 = (CF/CA-x*(1.-x))*t_channel_mD2->get_mD2(T)/2.;
  double Jacobian = 2*kT2/sintheta;
  double dR_dxdy = alpha_s(kT2, T)/(2.*M_PI)
                     * ( P_g2qqbar(x) + x*(1.-x)*m2*m2/(kT2+mg2+m2*m2))
                     / std::pow(kT2+mg2, 2)
                     * Jacobian * (CF/CA-x*(1.-x));
  return dR_dxdy;
}
