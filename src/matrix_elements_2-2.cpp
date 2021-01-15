#include <cmath>
#include <iostream>
#include "predefine.h"
#include "matrix_elements.h"
#include "lorentz.h"
#include "simpleLogger.h"

// parameters:
// s, tcut, T
// masses

//    g(h)1 + g(s)2 --only-t--> g3 + g4
double M2_gg2gg(const double t, void * params){
    double * p = static_cast<double*>(params);
    double s = p[0], tcut = p[1], Temp = p[2];
    //double m1=p[3], m2=p[4], m3=p[5], m4=p[6];
    if (t>tcut) return 0.;
    double Q2s = s, Q2t = t, Q2u = - s - t;
    double At = alpha_s(Q2t, Temp);
    double result = c72pi2*At*At*(-Q2s*Q2u/Q2t/Q2t);
    // There should also be a u channel divergence term, but it is 
    // accounted by taking the full space-space without identical
    // boson constrain into account.
    if (result < 0.) return 0.;
    return result;
}

double dX_gg2gg(const double t, void * params){
    double * p = static_cast<double*>(params);
    double s = p[0];
    return M2_gg2gg(t, params)/c16pi/s/s;
}

///    g(h)1 + q(s)2 --only 1/t^2 term--> g3 + q4
double M2_gq2gq(const double t, void * params){
    double * p = static_cast<double*>(params);
    double s = p[0], tcut = p[1], Temp = p[2];
    //double m1=p[3], m2=p[4], m3=p[5], m4=p[6];
    if (t>tcut) return 0;
    double Q2s = s, Q2t = t, Q2u = - s - t;
    double At = alpha_s(Q2t, Temp);
    double result = At*At*c16pi2*(Q2s*Q2s+Q2u*Q2u)/(Q2t*Q2t);
    if (result < 0.) return 0.;
    return result;
}
double dX_gq2gq(const double t, void * params){
    double * p = static_cast<double*>(params);
    double s = p[0];
    return M2_gq2gq(t, params)/c16pi/s/s;
}


///    q(h)1 + q(s)2 --only-t--> q3 + q4, m1=m3, m2=m4=0
double M2_qq2qq(const double t, void * params){
    double * p = static_cast<double*>(params);
    double s = p[0], tcut = p[1], Temp = p[2];
    double m1sq = std::pow(p[3], 2);
    if (t>tcut) return 0;
    double Q2s = s - m1sq, Q2t = t, Q2u = m1sq - s - t;
    double At = alpha_s(Q2t, Temp);
    double result = c64d9pi2*At*At*(Q2u*Q2u + Q2s*Q2s + 2.*m1sq*t)/Q2t/Q2t;
    if (result < 0.) return 0.;
    return result;
}
double dX_qq2qq(const double t, void * params){
    double * p = static_cast<double*>(params);
    double s = p[0], m1sq = std::pow(p[3], 2);
    return M2_qq2qq(t, params)/c16pi/std::pow(s-m1sq, 2);
}


///    q(h)1 + g(s)2 --only-t--> q3 + g4, m1=m3, m2=m4=0
double M2_qg2qg(const double t, void * params) {
    double * p = static_cast<double*>(params);
    double s = p[0], tcut = p[1], Temp = p[2];
    double m1sq = std::pow(p[3], 2);
    //double m2=p[4], m3=p[5], m4=p[6];
    if (t>tcut) return 0;
    double Q2s = s - m1sq, Q2t = t, Q2u = m1sq - s - t;
    double At = alpha_s(Q2t, Temp);
    double result = At*At*c16pi2*(Q2s*Q2s+Q2u*Q2u)/(Q2t*Q2t);
    if (result < 0.) return 0.;
    else return result;
}
double dX_qg2qg(const double t, void * params){
    double * p = static_cast<double*>(params);
    double s = p[0], m1sq = std::pow(p[3], 2);
    return M2_qg2qg(t, params)/c16pi/std::pow(s-m1sq, 2);
}

///    q(h)1 + g(s)2 --s+u+t--> q3 + g4, m1=m3, m2=m4=0
double M2_qg2qg_stu(const double t, void * params){
    double * p = static_cast<double*>(params);
    double s = p[0], tcut = p[1], Temp = p[2];
    double m1sq = std::pow(p[3], 2);
    //double m2=p[4], m3=p[5], m4=p[6];
    double Q2s = s - m1sq, Q2t = t, Q2u = m1sq - s - t;
    if (Q2t>tcut) return 0.;
    if (Q2u>tcut) return 0.;
    double At = alpha_s(Q2t, Temp),
           Au = alpha_s(Q2u, Temp),
           As = alpha_s(Q2s, Temp);

    double result = 2.*At*At * Q2s*(-Q2u)/Q2t/Q2t // t*t
    + c4d9*As*As*( Q2s*(-Q2u) + 2.*m1sq*(Q2s + 2.*m1sq) )/std::pow(Q2s, 2) // s*s
    + c4d9*Au*Au*( Q2s*(-Q2u) + 2.*m1sq*(Q2u + 2.*m1sq) )/std::pow(Q2u, 2) // u*u
    + c1d9*As*Au*m1sq*(4.*m1sq - t)/Q2s/(-Q2u)     // s*u
    + At*As*( Q2s*(-Q2u) + m1sq*(Q2s - Q2u) )/Q2t/Q2s    // t*s
    - At*Au*( Q2s*(-Q2u) - m1sq*(Q2s - Q2u) )/Q2t/(-Q2u); // t*u
    if (result < 0.) return 0.;
    return result*c16pi2;
}
double dX_qg2qg_stu(const double t, void * params){
    double * p = static_cast<double*>(params);
    double s = p[0], m1sq = std::pow(p[3], 2);
    return M2_qg2qg_stu(t, params)/c16pi/std::pow(s-m1sq, 2);
}

// g(h)1 + g(s)2 --s+t+u--> q3 + qbar4, m1=m2=0, m3=m4
double M2_gg2qqbar(const double t, void * params){
    double * p = static_cast<double*>(params);
    double s = p[0], tcut = p[1], Temp = p[2];
    double m3sq = std::pow(p[5], 2);
    //double m1 =p[3], m2=p[4], m4=p[6];
    double mg2 = t_channel_mD2->get_mD2(Temp)/2.;
    double S = s+mg2;
    double T = t - m3sq-mg2;
    double U = - S - T-mg2;
    // define coupling constant for each channel
    double As = alpha_s(S, Temp);
    double At = alpha_s(T, Temp);
    double Au = alpha_s(U, Temp);
    return std::pow(M_PI,2) * (
     12.*U*T/S/S*As*As 
   + 8./3.*(U/T+T/U)*At*Au
   - 16./3.*m3sq*((2.*m3sq+T)/T/T*At*At + (2*m3sq+U)/U/U*Au*Au)
   + 6.*(T+U)/S*As*As
   + 6.*m3sq/S*std::pow(T-U,2)/T/U*At*Au 
   - 2./3.*m3sq*(S-4.*m3sq)/T/U*At*Au
    );
}
double dX_gg2qqbar(const double t, void * params){
    double * p = static_cast<double*>(params);
    double s = p[0];
    return M2_gg2qqbar(t, params)/c16pi/std::pow(s, 2);
}


// q(h)1 + qbar(s)2 --s+t+u--> q3 + qbar4, m1=m2=0, m3=m4, different flavor
double M2_qqbar2qqbar_diff(const double t, void * params){
    double * p = static_cast<double*>(params);
    double s = p[0], tcut = p[1], Temp = p[2];
    double m3sq = std::pow(p[5], 2);
    //double m1 =p[3], m2=p[4], m4=p[6];
    double mg2 = t_channel_mD2->get_mD2(Temp)/2.;
    double S = s+mg2;
    double T = t - m3sq-mg2;
    double U = - S - T-mg2;

    // define coupling constant for each channel
    double As = alpha_s(S, Temp);
    return 64./9.*std::pow(M_PI*As,2)*(T*T+U*U+2.*m3sq*S)/std::pow(S,2);
}

double dX_qqbar2qqbar_diff(const double t, void * params){
    double * p = static_cast<double*>(params);
    double s = p[0];
    return M2_qqbar2qqbar_diff(t, params)/c16pi/std::pow(s, 2);
}





