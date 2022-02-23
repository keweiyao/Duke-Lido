#include <cmath>
#include <iostream>
#include "predefine.h"
#include "matrix_elements.h"
#include "lorentz.h"
#include "simpleLogger.h"


void unpack_three_body_phase_space(
    const double * x_, void * params_,
    double & kx, double & ky, double & yk,
    double & qx, double & qy,
    double & xbar, double & t,
    double & Jacobian, bool & status
    ){
    double * p = static_cast<double*>(params_);
    double s = p[0], tcut = p[1], Temp = p[2];
    double mA=p[3], mB=p[4], 
           m1=p[5], m2=p[6], m3=p[7];
    double sqrts = std::sqrt(s);
    double pT3_sq = x_[0],
           y3 = x_[1],
           costheta12 = x_[2], // in 1+2 com frame
           phi12 = x_[3];      // in 1+2 com frame
    double pT3 = std::sqrt(pT3_sq);
    double p3 = pT3*std::cosh(y3);
    double E3 = std::sqrt(p3*p3+m3*m3);
    double s12 = s+m3*m3-2.*sqrts*E3;
    bool out_of_bounds = std::fabs(costheta12)>=1. 
                      || phi12<=0. || phi12>=2.*M_PI 
                      || s12 >= s || s12 < std::pow(m1+m2, 2)
                      || s<=std::pow(mA+mB,2) || s<=std::pow(m1+m2+m3,2);
    if (out_of_bounds) { status = false; return;}
    else status = true;
    double EA = (s+mA*mA-mB*mB)/2./sqrts;
    double pA = std::sqrt(EA*EA-mA*mA);
    double EB = (s-mA*mA+mB*mB)/2./sqrts;
    double pB = std::sqrt(EB*EB-mB*mB);
    // construct k = p3
    fourvec kmu{E3, pT3, 0, pT3*std::sinh(y3)};
    // construct p1, p1 in (1+2)-CoM frame
    double cosphi12 = std::cos(phi12), sinphi12 = std::sin(phi12),
           sintheta12 = std::sqrt(1.-std::pow(costheta12,2));
    kx = kmu.x();
    ky = kmu.y();
    yk = y3;
    double V0 = sqrts - kmu.t();
    double v12[3] = {-kmu.x()/V0, -kmu.y()/V0, -kmu.z()/V0};
    double p1_12 = std::sqrt( (s12-std::pow(m1+m2,2)) 
                            * (s12-std::pow(m1-m2,2))
                            / 4. / s12 );
    Jacobian = p3/E3/c256pi4*p1_12/std::sqrt(s12);
     if (y3>=0){
        fourvec p2mu{std::sqrt(p1_12*p1_12+mB*mB), 
         -p1_12*sintheta12*cosphi12,
         -p1_12*sintheta12*sinphi12, 
         -p1_12*costheta12};   
        double m11 = m1+m2+m3-mB-m3;
        fourvec p1mu{std::sqrt(p1_12*p1_12+m11*m11), 
         p1_12*sintheta12*cosphi12,
         p1_12*sintheta12*sinphi12, 
         p1_12*costheta12};   
        p2mu = p2mu.boost_back(v12[0], v12[1], v12[2]);
        p1mu = p1mu.boost_back(v12[0], v12[1], v12[2]);
        qx = -p2mu.x();
        qy = -p2mu.y();
        t = 2.*mB*mB-2.*(EB*p2mu.t()+pB*p2mu.z());
        xbar = dot(kmu, p2mu)/( dot(kmu, p2mu) + dot(p1mu, p2mu) );
    }else{
        fourvec p1mu{std::sqrt(p1_12*p1_12+mA*mA), 
         p1_12*sintheta12*cosphi12,
         p1_12*sintheta12*sinphi12, 
         p1_12*costheta12};
        double m21 = m1+m2+m3-mA-m3;
        fourvec p2mu{std::sqrt(p1_12*p1_12+m21*m21), 
        -p1_12*sintheta12*cosphi12,
        -p1_12*sintheta12*sinphi12, 
        -p1_12*costheta12};   
        p1mu = p1mu.boost_back(v12[0], v12[1], v12[2]);
        p2mu = p2mu.boost_back(v12[0], v12[1], v12[2]);
        qx = -p1mu.x();
        qy = -p1mu.y();
        t = 2.*mA*mA-2.*(EA*p1mu.t()-pA*p1mu.z());
        xbar = dot(kmu, p1mu)/( dot(kmu, p1mu) + dot(p2mu, p1mu) );
    }   
}




void unpack_three_body_phase_space_same(
    const double * x_, void * params_,
    double & kx, double & ky, double & yk,
    double & qx, double & qy,
    double & xbar, double & t,
    double & Jacobian, bool & status
    ){
    double * p = static_cast<double*>(params_);
    double s = p[0], tcut = p[1], Temp = p[2];
    double mA=p[3], mB=p[4], 
           m1=p[5], m2=p[6], m3=p[7];
    double sqrts = std::sqrt(s);
    double pT3_sq = x_[0],
           y3 = x_[1],
           costheta12 = x_[2], // in 1+2 com frame
           phi12 = x_[3];      // in 1+2 com frame
    double pT3 = std::sqrt(pT3_sq);
    double p3 = pT3*std::cosh(y3);
    double E3 = std::sqrt(p3*p3+m3*m3);
    double s12 = s+m3*m3-2.*sqrts*E3;
    bool out_of_bounds = std::fabs(costheta12)>=1. 
                      || phi12<=0. || phi12>=2.*M_PI 
                      || s12 >= s || s12 < std::pow(m1+m2, 2)
                      || s<=std::pow(mA+mB,2) || s<=std::pow(m1+m2+m3,2);
    if (out_of_bounds) { status = false; return;}
    else status = true;
    double EA = (s+mA*mA-mB*mB)/2./sqrts;
    double pA = std::sqrt(EA*EA-mA*mA);
    double EB = (s-mA*mA+mB*mB)/2./sqrts;
    double pB = std::sqrt(EB*EB-mB*mB);
    // construct k = p3
    fourvec kmu{E3, pT3, 0, pT3*std::sinh(y3)};
    // construct p1, p1 in (1+2)-CoM frame
    double cosphi12 = std::cos(phi12), sinphi12 = std::sin(phi12),
           sintheta12 = std::sqrt(1.-std::pow(costheta12,2));
    kx = kmu.x();
    ky = kmu.y();
    yk = y3;
    double V0 = sqrts - kmu.t();
    double v12[3] = {-kmu.x()/V0, -kmu.y()/V0, -kmu.z()/V0};
    double p1_12 = std::sqrt( (s12-std::pow(m1+m2,2)) 
                            * (s12-std::pow(m1-m2,2))
                            / 4. / s12 );
    Jacobian = p3/E3/c256pi4*p1_12/std::sqrt(s12);
fourvec p2mu{std::sqrt(p1_12*p1_12+mB*mB), 
         -p1_12*sintheta12*cosphi12,
         -p1_12*sintheta12*sinphi12, 
         -p1_12*costheta12};   
        double m11 = m1+m2+m3-mB-m3;
        fourvec p1mu{std::sqrt(p1_12*p1_12+m11*m11), 
         p1_12*sintheta12*cosphi12,
         p1_12*sintheta12*sinphi12, 
         p1_12*costheta12};   
        p2mu = p2mu.boost_back(v12[0], v12[1], v12[2]);
        p1mu = p1mu.boost_back(v12[0], v12[1], v12[2]);
        qx = -p2mu.x();
        qy = -p2mu.y();
        t = 2.*mB*mB-2.*(EB*p2mu.t()+pB*p2mu.z());
        xbar = dot(kmu, p2mu)/( dot(kmu, p2mu) + dot(p1mu, p2mu) );

}
// q[h](A) + q[s](B) --> q(1)+q(2)+g(3)
double M2_qq2qqg(const double * x_, void * params_){
    double * params = static_cast<double*>(params_);
    double Temp = params[2];
    double kx, ky, qx, qy, xbar, t, Jacobian, yk;
    bool status = false;
    unpack_three_body_phase_space(
      x_, params, 
      kx, ky, yk, 
      qx, qy,
      xbar, t,
      Jacobian, status
    );
    if (status==false) return 0.;
    double one_minus_xbar = 1.-xbar;
    double mg2 = t_channel_mD2->get_mD2(Temp) / 2.;
    double mth2 = (CF/CA*xbar*xbar+1.-xbar)*mg2;
    double kt2 = kx*kx+ky*ky;
    double Msq_22 = M2_qq2qq(t, params);
    double MA[2], MB[2], MC[2], A[2], B[2];
    MA[0] = kx - qx;
    MA[1] = ky - qy;
    MB[0] = kx - xbar*qx;
    MB[1] = ky - xbar*qy;
    MC[0] = kx;
    MC[1] = ky;  
    double DA = MA[0]*MA[0] + MA[1]*MA[1] + mth2;
    double DB = MB[0]*MB[0] + MB[1]*MB[1] + mth2;
    double DC = MC[0]*MC[0] + MC[1]*MC[1] + mth2;
    A[0] = MA[0]/DA-MB[0]/DB;
    A[1] = MA[1]/DA-MB[1]/DB;
    B[0] = MA[0]/DA-MC[0]/DC;
    B[1] = MA[1]/DA-MC[1]/DC;  
    double A2 = std::pow(A[0],2) + std::pow(A[1],2);
    double B2 = std::pow(B[0],2) + std::pow(B[1],2);
    double AB = A[0]*B[0] + A[1]*B[1];
    double Msq_12 = c8pi * alpha_s(kt2, Temp) 
            * xbar * one_minus_xbar
            * P_q2qg(xbar)
            * (CF*A2 + CF*B2 - (2.*CF-CA)*AB) / CF;

    // 2->3 = 2->2 * 1->2
    return Msq_22*Msq_12*Jacobian;
}
double dX_qq2qqg(const double * x_, void * params){
    double * p = static_cast<double*>(params);
    double s = p[0];
    double mA=p[3], mB=p[4];
    double Flux = 2*std::sqrt((s-std::pow(mA-mB,2))*(s-std::pow(mA+mB,2)));
    return M2_qq2qqg(x_, params)/Flux;
}


// q[h](A) + g[s](B) --> q(1)+g(2)+g(3)
double M2_qg2qgg(const double * x_, void * params_){
    double * params = static_cast<double*>(params_);
    double Temp = params[2];
    double kx, ky, qx, qy, xbar, t, Jacobian, yk;
    bool status = false;
    unpack_three_body_phase_space(
      x_, params, 
      kx, ky, yk, 
      qx, qy,
      xbar, t,
      Jacobian, status
    );
    if (status==false) return 0.;
    double one_minus_xbar = 1.-xbar;
    double kt2 = kx*kx+ky*ky;
    double MA[2], MB[2], MC[2], A[2], B[2];
    double DA, DB, DC, A2, B2, AB, Msq_22, Msq_12;
    if (yk>0){
    double mg2 = t_channel_mD2->get_mD2(Temp) / 2.;
    double mth2 = (CF/CA*xbar*xbar+1.-xbar)*mg2;
        Msq_22 = M2_qg2qg(t, params);
        MA[0] = kx - qx;
        MA[1] = ky - qy;
        MB[0] = kx - xbar*qx;
        MB[1] = ky - xbar*qy;
        MC[0] = kx;
        MC[1] = ky;
        DA = MA[0]*MA[0] + MA[1]*MA[1] + mth2;
        DB = MB[0]*MB[0] + MB[1]*MB[1] + mth2;
        DC = MC[0]*MC[0] + MC[1]*MC[1] + mth2;
        A[0] = MA[0]/DA-MB[0]/DB;
        A[1] = MA[1]/DA-MB[1]/DB;
        B[0] = MA[0]/DA-MC[0]/DC;
        B[1] = MA[1]/DA-MC[1]/DC;   
        A2 = std::pow(A[0],2) + std::pow(A[1],2);
        B2 = std::pow(B[0],2) + std::pow(B[1],2);
        AB = A[0]*B[0] + A[1]*B[1];     
        Msq_12 = c8pi* alpha_s(kt2, Temp)
            * xbar * one_minus_xbar
            * P_q2qg(xbar)
            * (CF*A2 + CF*B2 - (2.*CF-CA)*AB) /CF;
    }
    else{
        if (xbar>.5) return 0.;
    double mg2 = t_channel_mD2->get_mD2(Temp) / 2.;
    double mth2 = (xbar*xbar+1.-xbar)*mg2;
        Msq_22 = M2_qg2qg(t, params);
        MA[0] = kx - xbar*qx;
        MA[1] = ky - xbar*qy;
        MB[0] = kx - qx;
        MB[1] = ky - qy;
        MC[0] = kx;
        MC[1] = ky;
        DA = MA[0]*MA[0] + MA[1]*MA[1] + mth2;
        DB = MB[0]*MB[0] + MB[1]*MB[1] + mth2;
        DC = MC[0]*MC[0] + MC[1]*MC[1] + mth2;
        A[0] = MA[0]/DA-MB[0]/DB;
        A[1] = MA[1]/DA-MB[1]/DB;
        B[0] = MA[0]/DA-MC[0]/DC;
        B[1] = MA[1]/DA-MC[1]/DC;
        A2 = std::pow(A[0],2) + std::pow(A[1],2);
        B2 = std::pow(B[0],2) + std::pow(B[1],2);
        AB = A[0]*B[0] + A[1]*B[1];
        Msq_12 = c8pi* alpha_s(kt2, Temp)
            * xbar * one_minus_xbar
            * P_g2gg(xbar)
            * (A2 + B2 - AB);
    }
    return Msq_12*Msq_22*Jacobian;
}    
double dX_qg2qgg(const double * x_, void * params){
    double * p = static_cast<double*>(params);
    double s = p[0];
    double mA=p[3], mB=p[4];
    double Flux = 2*std::sqrt((s-std::pow(mA-mB,2))*(s-std::pow(mA+mB,2)));
    return M2_qg2qgg(x_, params)/Flux;
}

// q[h](A) + g[s](B) --> q(1)+q(2)+qbar(3)
double M2_qg2qqqbar(const double * x_, void * params_){
    double * params = static_cast<double*>(params_);
    double Temp = params[2], m3 = params[7];
    double kx, ky, qx, qy, xbar, t, Jacobian, yk;
    bool status = false;
    unpack_three_body_phase_space(
      x_, params, 
      kx, ky, yk, 
      qx, qy,
      xbar, t,
      Jacobian, status
    );
    if (status==false) return 0.;
    double one_minus_xbar = 1.-xbar;
    double mg2 = t_channel_mD2->get_mD2(Temp) / 2.;
    double mth2 = (CF/CA-xbar*(1-xbar))*mg2 + m3*m3;
    double kt2 = kx*kx+ky*ky;
    double MA[2], MB[2], MC[2], A[2], B[2];
    double DA, DB, DC, A2, B2, AB, Msq_22, Msq_12;
    Msq_22 = M2_qg2qg(t, params);
    if (yk>0){
        return 0.;
    }
    else{
        MA[0] = kx - xbar*qx;
        MA[1] = ky - xbar*qy;
        MB[0] = kx - qx;
        MB[1] = ky - qy;
        MC[0] = kx;
        MC[1] = ky;
        DA = MA[0]*MA[0] + MA[1]*MA[1] + mth2;
        DB = MB[0]*MB[0] + MB[1]*MB[1] + mth2;
        DC = MC[0]*MC[0] + MC[1]*MC[1] + mth2;
        A[0] = MA[0]/DA-MB[0]/DB;
        A[1] = MA[1]/DA-MB[1]/DB;
        B[0] = MA[0]/DA-MC[0]/DC;
        B[1] = MA[1]/DA-MC[1]/DC;
        A2 = std::pow(A[0],2) + std::pow(A[1],2);
        B2 = std::pow(B[0],2) + std::pow(B[1],2);
        AB = A[0]*B[0] + A[1]*B[1];
        double A2M = std::pow(1/DA-1/DB,2);
        double B2M = std::pow(1/DA-1/DC,2);
        double ABM = (1/DA-1/DB)*(1/DA-1/DC);
        Msq_12 = c8pi* alpha_s(kt2, Temp) 
            * one_minus_xbar * xbar
            * ( P_g2qqbar(xbar)*( CF*A2 + CF*B2 - (2*CF-CA)*AB ) // anti-spin
              + m3*m3/2. * ( CF*A2M + CF*B2M - (2*CF-CA)*ABM ) // para-spin
            ) /CA;

    }
    return Msq_12*Msq_22*Jacobian;
}    
double dX_qg2qqqbar(const double * x_, void * params){
    double * p = static_cast<double*>(params);
    double s = p[0];
    double mA=p[3], mB=p[4];
    double Flux = 2*std::sqrt((s-std::pow(mA-mB,2))*(s-std::pow(mA+mB,2)));
    return M2_qg2qqqbar(x_, params)/Flux;
}


// g[h](A) + q[s](B) --> g(1)+q(2)+g(3)
double M2_gq2gqg(const double * x_, void * params_){
    double * params = static_cast<double*>(params_);
    double Temp = params[2];
    double kx, ky, qx, qy, xbar, t, Jacobian, yk;
    bool status = false;
    unpack_three_body_phase_space(
      x_, params, 
      kx, ky, yk, 
      qx, qy,
      xbar, t,
      Jacobian, status
    );
    if (status==false) return 0.;
    double one_minus_xbar = 1.-xbar;
    double kt2 = kx*kx+ky*ky;
    double MA[2], MB[2], MC[2], A[2], B[2];
    double DA, DB, DC, A2, B2, AB, Msq_22, Msq_12;
    Msq_22 = M2_gq2gq(t, params);
    if (yk>0){
        if (xbar>.5) return 0.;
    double mg2 = t_channel_mD2->get_mD2(Temp) / 2.;
    double mth2 = (1.-xbar+xbar*xbar)*mg2;
        MA[0] = kx - xbar*qx;
        MA[1] = ky - xbar*qy;
        MB[0] = kx - qx;
        MB[1] = ky - qy;
        MC[0] = kx;
        MC[1] = ky;
        DA = MA[0]*MA[0] + MA[1]*MA[1] + mth2;
        DB = MB[0]*MB[0] + MB[1]*MB[1] + mth2;
        DC = MC[0]*MC[0] + MC[1]*MC[1] + mth2;
        A[0] = MA[0]/DA-MB[0]/DB;
        A[1] = MA[1]/DA-MB[1]/DB;
        B[0] = MA[0]/DA-MC[0]/DC;
        B[1] = MA[1]/DA-MC[1]/DC;
        A2 = std::pow(A[0],2) + std::pow(A[1],2);
        B2 = std::pow(B[0],2) + std::pow(B[1],2);
        AB = A[0]*B[0] + A[1]*B[1];
        Msq_12 = c8pi* alpha_s(kt2, Temp)
            * xbar * one_minus_xbar
            *  P_g2gg(xbar)
            * (A2 + B2 - AB);
    }
    else{
    double mg2 = t_channel_mD2->get_mD2(Temp) / 2.;
    double mth2 = (CF/CA*xbar*xbar+1-xbar)*mg2;
        MA[0] = kx - qx;
        MA[1] = ky - qy;
        MB[0] = kx - xbar*qx;
        MB[1] = ky - xbar*qy;
        MC[0] = kx;
        MC[1] = ky;
        DA = MA[0]*MA[0] + MA[1]*MA[1] + mth2;
        DB = MB[0]*MB[0] + MB[1]*MB[1] + mth2;
        DC = MC[0]*MC[0] + MC[1]*MC[1] + mth2;
        A[0] = MA[0]/DA-MB[0]/DB;
        A[1] = MA[1]/DA-MB[1]/DB;
        B[0] = MA[0]/DA-MC[0]/DC;
        B[1] = MA[1]/DA-MC[1]/DC;   
        A2 = std::pow(A[0],2) + std::pow(A[1],2);
        B2 = std::pow(B[0],2) + std::pow(B[1],2);
        AB = A[0]*B[0] + A[1]*B[1];     
        Msq_12 = c8pi* alpha_s(kt2, Temp) 
            * xbar * one_minus_xbar
            * P_q2qg(xbar)
            * (CF*A2 + CF*B2 - (2.*CF-CA)*AB)/6;
    }
    return Msq_12*Msq_22*Jacobian;
}    
double dX_gq2gqg(const double * x_, void * params){
    double * p = static_cast<double*>(params);
    double s = p[0];
    double mA=p[3], mB=p[4];
    double Flux = 2*std::sqrt((s-std::pow(mA-mB,2))*(s-std::pow(mA+mB,2)));
    return M2_gq2gqg(x_, params)/Flux;
}



// g[h](A) + g[s](B) --> g(1)+g(2)+g(3)
double M2_gg2ggg(const double * x_, void * params_){
    double * params = static_cast<double*>(params_);
    double Temp = params[2];
    double kx, ky, qx, qy, xbar, t, Jacobian, yk;
    bool status = false;
    unpack_three_body_phase_space(
      x_, params, 
      kx, ky, yk, 
      qx, qy,
      xbar, t,
      Jacobian, status
    );
    if (status==false) return 0.;
    if (xbar>.5) return 0.;

    double one_minus_xbar = 1.-xbar;
    double mg2 = t_channel_mD2->get_mD2(Temp) / 2.;
    double mth2 = (1.-xbar+xbar*xbar)*mg2;
    double kt2 = kx*kx+ky*ky;
    double MA[2], MB[2], MC[2], A[2], B[2];
    double DA, DB, DC, A2, B2, AB, Msq_22, Msq_12;
    
    MA[0] = kx - xbar*qx;
    MA[1] = ky - xbar*qy;
    MB[0] = kx - qx;
    MB[1] = ky - qy;
    MC[0] = kx;
    MC[1] = ky;
    DA = MA[0]*MA[0] + MA[1]*MA[1] + mth2;
    DB = MB[0]*MB[0] + MB[1]*MB[1] + mth2;
    DC = MC[0]*MC[0] + MC[1]*MC[1] + mth2;
    A[0] = MA[0]/DA-MB[0]/DB;
    A[1] = MA[1]/DA-MB[1]/DB;
    B[0] = MA[0]/DA-MC[0]/DC;
    B[1] = MA[1]/DA-MC[1]/DC;
    A2 = std::pow(A[0],2) + std::pow(A[1],2);
    B2 = std::pow(B[0],2) + std::pow(B[1],2);
    AB = A[0]*B[0] + A[1]*B[1];
    Msq_22 = M2_gg2gg(t, params);
    Msq_12 = c8pi* alpha_s(kt2, Temp)
        * xbar * one_minus_xbar
        * P_g2gg(xbar)
        * (A2 + B2 - AB);
    return Msq_12*Msq_22*Jacobian;
}    
double dX_gg2ggg(const double * x_, void * params){
    double * p = static_cast<double*>(params);
    double s = p[0];
    double mA=p[3], mB=p[4];
    double Flux = 2*std::sqrt((s-std::pow(mA-mB,2))*(s-std::pow(mA+mB,2)));
    return M2_gg2ggg(x_, params)/Flux;
}

// g[h](A) + q[s](B) --> q(1)+q(2)+qbar(3)
double M2_gq2qqqbar(const double * x_, void * params_){
    double * params = static_cast<double*>(params_);
    double Temp = params[2], m3 = params[7];
    double kx, ky, qx, qy, xbar, t, Jacobian, yk;
    bool status = false;
    unpack_three_body_phase_space_same(
      x_, params, 
      kx, ky, yk, 
      qx, qy,
      xbar, t,
      Jacobian, status
    );
    if (status==false) return 0.;
    double one_minus_xbar = 1.-xbar;
    double mg2 = t_channel_mD2->get_mD2(Temp) / 2.;
    double kt2 = kx*kx+ky*ky;
    double MA[2], MB[2], MC[2], A[2], B[2];
    double DA, DB, DC, A2, B2, AB, Msq_22, Msq_12;
    if (yk<0.) return 0.;
    MA[0] = kx - xbar*qx;
    MA[1] = ky - xbar*qy;
    MB[0] = kx - qx;
    MB[1] = ky - qy;
    MC[0] = kx;
    MC[1] = ky;
    double mth2 = (CF/CA-xbar+xbar*xbar)*mg2 + m3*m3;
    DA = MA[0]*MA[0] + MA[1]*MA[1] + mth2;
    DB = MB[0]*MB[0] + MB[1]*MB[1] + mth2;
    DC = MC[0]*MC[0] + MC[1]*MC[1] + mth2;
    A[0] = MA[0]/DA-MB[0]/DB;
    A[1] = MA[1]/DA-MB[1]/DB;
    B[0] = MA[0]/DA-MC[0]/DC;
    B[1] = MA[1]/DA-MC[1]/DC;
    A2 = std::pow(A[0],2) + std::pow(A[1],2);
    B2 = std::pow(B[0],2) + std::pow(B[1],2);
    AB = A[0]*B[0] + A[1]*B[1];
    double A2M = std::pow(1/DA-1/DB,2);
    double B2M = std::pow(1/DA-1/DC,2);
    double ABM = (1/DA-1/DB)*(1/DA-1/DC);
    Msq_22 = M2_gq2gq(t, params);
    Msq_12 = c8pi* alpha_s(kt2, Temp) 
            * one_minus_xbar * xbar
            * ( P_g2qqbar(xbar)*( CF*A2 + CF*B2 - (2*CF-CA)*AB ) // anti-spin
              + m3*m3/2. * ( CF*A2M + CF*B2M - (2*CF-CA)*ABM ) // para-spin
            ) /CA;

    return Msq_12*Msq_22*Jacobian;
}    
double dX_gq2qqqbar(const double * x_, void * params){
    double * p = static_cast<double*>(params);
    double s = p[0];
    double mA=p[3], mB=p[4];
    double Flux = 2*std::sqrt((s-std::pow(mA-mB,2))*(s-std::pow(mA+mB,2)));
    return M2_gq2qqqbar(x_, params)/Flux;
}

// g[h](A) + g[s](B) --> q(1)+g(2)+qbar(3)
double M2_gg2qgqbar(const double * x_, void * params_){
    double * params = static_cast<double*>(params_);
    double Temp = params[2], m3 = params[7];
    double kx, ky, qx, qy, xbar, t, Jacobian, yk;
    bool status = false;
    unpack_three_body_phase_space(
      x_, params, 
      kx, ky, yk, 
      qx, qy,
      xbar, t,
      Jacobian, status
    );
    if (status==false) return 0.;
    if (yk<0) return 0.;
    double one_minus_xbar = 1.-xbar;
    double mg2 = t_channel_mD2->get_mD2(Temp) / 2.;
    double kt2 = kx*kx+ky*ky;
    double MA[2], MB[2], MC[2], A[2], B[2];
    double DA, DB, DC, A2, B2, AB, Msq_22, Msq_12;
    MA[0] = kx - xbar*qx;
    MA[1] = ky - xbar*qy;
    MB[0] = kx - qx;
    MB[1] = ky - qy;
    MC[0] = kx;
    MC[1] = ky;
    double mth2 = (CF/CA-xbar+xbar*xbar)*mg2 + m3*m3;
    DA = MA[0]*MA[0] + MA[1]*MA[1] + mth2;
    DB = MB[0]*MB[0] + MB[1]*MB[1] + mth2;
    DC = MC[0]*MC[0] + MC[1]*MC[1] + mth2;
    A[0] = MA[0]/DA-MB[0]/DB;
    A[1] = MA[1]/DA-MB[1]/DB;
    B[0] = MA[0]/DA-MC[0]/DC;
    B[1] = MA[1]/DA-MC[1]/DC;
    A2 = std::pow(A[0],2) + std::pow(A[1],2);
    B2 = std::pow(B[0],2) + std::pow(B[1],2);
    AB = A[0]*B[0] + A[1]*B[1];
    double A2M = std::pow(1/DA-1/DB,2);
    double B2M = std::pow(1/DA-1/DC,2);
    double ABM = (1/DA-1/DB)*(1/DA-1/DC);
    Msq_22 = M2_gg2gg(t, params);
    Msq_12 = c8pi* alpha_s(kt2, Temp) 
            * one_minus_xbar * xbar
            * ( P_g2qqbar(xbar)*( CF*A2 + CF*B2 - (2*CF-CA)*AB ) // anti-spin
              + m3*m3/2. * ( CF*A2M + CF*B2M - (2*CF-CA)*ABM ) // para-spin
            ) /CA;
    return Msq_12*Msq_22*Jacobian;
}    
double dX_gg2qgqbar(const double * x_, void * params){
    double * p = static_cast<double*>(params);
    double s = p[0];
    double mA=p[3], mB=p[4];
    double Flux = 2*std::sqrt((s-std::pow(mA-mB,2))*(s-std::pow(mA+mB,2)));
    return M2_gg2qgqbar(x_, params)/Flux;
}


// g[h](A) + g[s](B) --> q(1)+g(2)+qbar(3)
double M2_gg2gqqbar(const double * x_, void * params_){
    double * params = static_cast<double*>(params_);
    double Temp = params[2], m3 = params[7];
    double kx, ky, qx, qy, xbar, t, Jacobian, yk;
    bool status = false;
    unpack_three_body_phase_space(
      x_, params, 
      kx, ky, yk, 
      qx, qy,
      xbar, t,
      Jacobian, status
    );
    if (status==false) return 0.;
    if (yk>0) return 0.;
    double one_minus_xbar = 1.-xbar;
    double mg2 = t_channel_mD2->get_mD2(Temp) / 2.;
    double kt2 = kx*kx+ky*ky;
    double MA[2], MB[2], MC[2], A[2], B[2];
    double DA, DB, DC, A2, B2, AB, Msq_22, Msq_12;
    MA[0] = kx - xbar*qx;
    MA[1] = ky - xbar*qy;
    MB[0] = kx - qx;
    MB[1] = ky - qy;
    MC[0] = kx;
    MC[1] = ky;
    double mth2 = (CF/CA-xbar+xbar*xbar)*mg2 + m3*m3;
    DA = MA[0]*MA[0] + MA[1]*MA[1] + mth2;
    DB = MB[0]*MB[0] + MB[1]*MB[1] + mth2;
    DC = MC[0]*MC[0] + MC[1]*MC[1] + mth2;
    A[0] = MA[0]/DA-MB[0]/DB;
    A[1] = MA[1]/DA-MB[1]/DB;
    B[0] = MA[0]/DA-MC[0]/DC;
    B[1] = MA[1]/DA-MC[1]/DC;
    A2 = std::pow(A[0],2) + std::pow(A[1],2);
    B2 = std::pow(B[0],2) + std::pow(B[1],2);
    AB = A[0]*B[0] + A[1]*B[1];
    double A2M = std::pow(1/DA-1/DB,2);
    double B2M = std::pow(1/DA-1/DC,2);
    double ABM = (1/DA-1/DB)*(1/DA-1/DC);
    Msq_22 = M2_gg2gg(t, params);
    Msq_12 = c8pi* alpha_s(kt2, Temp) 
            * one_minus_xbar * xbar
            * ( P_g2qqbar(xbar)*( CF*A2 + CF*B2 - (2*CF-CA)*AB ) // anti-spin
	      + m3*m3/2. * ( CF*A2M + CF*B2M - (2*CF-CA)*ABM ) // para-spin
	    ) /CA;
    return Msq_12*Msq_22*Jacobian;
}    
double dX_gg2gqqbar(const double * x_, void * params){
    double * p = static_cast<double*>(params);
    double s = p[0];
    double mA=p[3], mB=p[4];
    double Flux = 2*std::sqrt((s-std::pow(mA-mB,2))*(s-std::pow(mA+mB,2)));
    return M2_gg2gqqbar(x_, params)/Flux;
}



