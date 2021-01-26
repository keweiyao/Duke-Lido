#include "Xsection.h"
#include "matrix_elements.h"
#include "integrator.h"
#include "minimizer.h"
#include "sampler.h"
#include "predefine.h"
#include <fstream>
#include "random.h"
#include "approx_functions.h"

template<>
Xsection<HS2PP, 2, double(*)(const double, void*)>::
    Xsection(std::string Name, std::string configfile,
            double(*f)(const double, void*)):
StochasticBase<2>(Name+"/xsection", configfile),
_f(f)
{
    // read configfile
    boost::property_tree::ptree config;
    std::ifstream input(configfile);
    read_xml(input, config);

    std::vector<std::string> strs;
    boost::split(strs, Name, boost::is_any_of("/"));
    std::string model_name = strs[0];
    std::string process_name = strs[1];
    auto tree = config.get_child(model_name+"."+process_name);
    _process_id = get_process_info(process_name, _IS_masses, _FS_masses, 
                                   _IS_types, _FS_types);

    // Set Approximate function for X and dX_max
    StochasticBase<2>::_ZeroMoment->SetApproximateFunction(approx_X22);
    StochasticBase<2>::_FunctionMax->SetApproximateFunction(approx_dX22_max);
}

template<>
Xsection<HS2PPP, 2, double(*)(const double *, void*)>::
    Xsection(std::string Name, std::string configfile,
            double(*f)(const double*, void*)):
StochasticBase<2>(Name+"/xsection", configfile),
_f(f)
{
    // read configfile
    boost::property_tree::ptree config;
    std::ifstream input(configfile);
    read_xml(input, config);

    std::vector<std::string> strs;
    boost::split(strs, Name, boost::is_any_of("/"));
    std::string model_name = strs[0];
    std::string process_name = strs[1];
    auto tree = config.get_child(model_name+"."+process_name);
    _process_id = get_process_info(process_name, _IS_masses, _FS_masses, 
                                   _IS_types, _FS_types);
    
}


/*****************************************************************/
/*************************sample dX/dPS **************************/
/*****************************************************************/
/*------------------Implementation for 2 -> 2--------------------*/
template<>
bool Xsection<HS2PP, 2, double(*)(const double, void*)>::
        sample(std::vector<double> parameters,
               int incoming_hard_pid,
               std::vector<fourvec> & FS,
               std::vector<int> & pids){
    double lnsqrts = parameters[0], temp = parameters[1];
    double sqrts = std::exp(lnsqrts);
    double s = std::pow(sqrts,2);
    double mA = _IS_masses[0], mB = _IS_masses[1];
    double m1 = _FS_masses[0], m2 = _FS_masses[1];
    double mD2 = t_channel_mD2->get_mD2(temp);
    double tcut = isPairProduction(_process_id)?
                -mD2/16.:-cut*mD2;
    if (s<=std::pow(std::sqrt(mA*mA-tcut/2.)+std::sqrt(mB*mB-tcut/2.), 2)
      ||s<=std::pow(std::sqrt(m1*m1-tcut/2.)+std::sqrt(m2*m2-tcut/2.), 2)) {
        LOG_INFO << "sample failed (below threshold) in 2->2 X";
	FS.clear();
        pids.clear();
        return false;  
    }
    double pAcm = std::sqrt(
                  (std::pow(sqrts-mA,2)-mB*mB)
                * (std::pow(sqrts+mA,2)-mB*mB)
                  /4./s);
    double p1cm = std::sqrt(
                  (std::pow(sqrts-m1,2)-m2*m2)
                * (std::pow(sqrts+m1,2)-m2*m2)
                  /4./s);

    double tm2 = std::pow((mA*mA-m1*m1-mB*mB+m2*m2)/2./sqrts, 2);
    double tmin0 = tm2 - std::pow(pAcm+p1cm, 2);
    double tmax0 = tm2 - std::pow(pAcm-p1cm, 2);
    double tmin = tmin0-tcut;
    double tmax = tmax0+tcut;
    double cosmax = 1. - (tmax0-tmax)/(2.*pAcm*p1cm);
    double cosmin = 1. - (tmax0-tmin)/(2.*pAcm*p1cm);
    double etamax = .5*std::log((1.+cosmax)/(1.-cosmax));
    double etamin = .5*std::log((1.+cosmin)/(1.-cosmin));
    // Define a lamda function, which transforms the integral
    // varaible for more efficient integration
    auto dXdeta = [s, temp, tmax, mA, mB, m1, m2, pAcm, p1cm, tmax0, this]
       (double eta) {
        double params[7] = {s, tmax, temp, mA, mB, m1, m2};
        double costheta = std::tanh(eta);
        double t = tmax0 - 2.*pAcm*p1cm*(1.-costheta);
        double Jacobian = 2.*pAcm*p1cm*(1.-costheta*costheta);
        return this->_f(t, params)*Jacobian;
    };
    double params[7] = {s, tmax, temp, mA, mB, m1, m2};
    double fmax = std::exp(StochasticBase<2>::GetFmax(parameters).s);
    bool status=true;
    double eta = sample_1d(dXdeta, {etamin, etamax}, fmax, status);
    if (status==true){
        double costheta = std::tanh(eta);
        double sintheta = std::sqrt(1.-costheta*costheta);
        double phi = Srandom::dist_phi(Srandom::gen);
        double cosphi = std::cos(phi), sinphi = std::sin(phi);
        double E1cm = (s+m1*m1-m2*m2)/2./sqrts;
        double E2cm = (s+m2*m2-m1*m1)/2./sqrts;

        FS.clear();
        FS.resize(2);
        FS[0] = {E1cm, p1cm*sintheta*cosphi, p1cm*sintheta*sinphi, p1cm*costheta};
        FS[1] = {E2cm, -p1cm*sintheta*cosphi, -p1cm*sintheta*sinphi, -p1cm*costheta};
        //LOG_INFO << "look: " << FS[0] << " " << FS[1] << ", " << cosmin<< " " << cosmax;
        assign_2to2_pid(_process_id, incoming_hard_pid, pids);
        return true;
    }else{
        LOG_INFO << "sample failed in 2->2 X";
	FS.clear();
        pids.clear();
        return false; 
    }
}



/*------------------Implementation for 2 -> 3--------------------*/
template<>
bool Xsection<HS2PPP, 2, double(*)(const double*, void*)>::
    sample(std::vector<double> parameters, 
           int incoming_hard_pid,
            std::vector<fourvec> & FS, 
            std::vector<int> & pids){
    double lnsqrts = parameters[0], temp = parameters[1];
    double mD2 = t_channel_mD2->get_mD2(temp);
    double tcut = is2to2_nlo(_process_id)?-mD2:-cut*mD2;
    double mA = _IS_masses[0], mB = _IS_masses[1];
    double m1 = _FS_masses[0], m2 = _FS_masses[1], m3 = _FS_masses[2];
    double sqrts = std::exp(lnsqrts);
    double s = std::pow(sqrts,2);
    if (s<std::pow(mA+mB,2)||s<std::pow(m1+m2+m3,2)) {
        LOG_INFO << "sample failed (below threshold) in 2->3 X";
	FS.clear();
        pids.clear();
        return false;
    }
    sqrts = std::sqrt(s);
    double pAcm = std::sqrt(
                  (std::pow(sqrts-mA,2)-mB*mB)
                * (std::pow(sqrts+mA,2)-mB*mB)
                  /4./s);

    //       pT3_sq 
    //       y3
    //       costheta12 // in 1+2 com frame
    //       phi12      // in 1+2 com frame
    // x0 = log(1+kt^2/mD2), 
    // x1 = y3/log(sqrts/kt)
    // x2 = log(1+2*(1-costheta12)*pAcm^2/mD2) ~ log(1+qt_cm12^2/mD2)
    // x3 = phi12
    // Jacobian = d(kT2, y3, c12, phi12)/d(x0, x1, x3, x3) = (kT2+mD2)*(c12+mD/2/pA/pA)*log(sqrts/kt)
    double mD2_2pA2 = mD2 / 2./pAcm/pAcm;
    auto dXdPS = [s, tcut, temp, mA, mB, m1, m2, m3, mD2, mD2_2pA2, pAcm, this](double * x_){
        double params[8] = {s, tcut, temp, mA, mB, m1, m2, m3};
        double kt2 = mD2*(std::exp(x_[0])-1);
        double ymax = std::acosh(pAcm/std::sqrt(kt2));
        double y3 = x_[1]*ymax, phi12 = x_[3]; 
        double costheta12 = 1.-mD2_2pA2*(std::exp(x_[2])-1);
        double x[4] = {kt2, y3, costheta12, phi12};
        double Jacobian = (1.-costheta12+ mD2_2pA2)*(kt2+mD2)*ymax;
        return this->_f(x, params)*Jacobian;
    };
    double xmin[4] = {0., -1.,  0., 0.};
    double xmax[4] = {std::log(1.+pAcm*pAcm/mD2), 1., std::log(1+2/mD2_2pA2), c2pi};
    double fmax = std::exp(StochasticBase<2>::GetFmax(parameters).s);
    bool status = true;
    auto res = sample_nd(dXdPS, 4, {{xmin[0], xmax[0]}, {xmin[1], xmax[1]},
                                    {xmin[2], xmax[2]}, {xmin[3], xmax[3]}},
                                    fmax, status);
    if (status==true) {
        // deconvolve parameter
        double kt2 = mD2*(std::exp(res[0])-1);
        double ymax = std::acosh(pAcm/std::sqrt(kt2));
        double yk = res[1]*ymax;
        double phi12 = res[3]; 
        double kt = std::sqrt(kt2);
        double costheta12 = 1.-mD2_2pA2*(std::exp(res[2])-1);
        assign_2to3_pid(_process_id, incoming_hard_pid, pids, (yk>=0));

        // reconstruct momentums
        double M1sq = std::pow(pid2mass(pids[0]),2);
        double M2sq = std::pow(pid2mass(pids[1]),2);
        double Mksq = std::pow(pid2mass(pids[2]),2);
        double Ek = std::sqrt(std::pow(kt*std::cosh(yk),2) + Mksq);
        fourvec kmu{Ek, kt, 0., kt*std::sinh(yk)}; 
        double s12 = s + m3*m3 - 2.*sqrts*Ek;
        double V0 = sqrts - kmu.t();
        double v12[3] = {-kmu.x()/V0, -kmu.y()/V0, -kmu.z()/V0};
        double p1_12cm = std::sqrt( (s12-std::pow(m1+m2,2)) 
                                  * (s12-std::pow(m1-m2,2))
                                  / 4. / s12 );
        // construct p1, p2 in (1+2)-com frame
        double cosphi12 = std::cos(phi12), sinphi12 = std::sin(phi12),
               sintheta12 = std::sqrt(1.-std::pow(costheta12,2));
        fourvec p1mu = fourvec{std::sqrt(p1_12cm*p1_12cm+M1sq), p1_12cm*sintheta12*cosphi12, p1_12cm*sintheta12*sinphi12, p1_12cm*costheta12};
        fourvec p2mu = fourvec{std::sqrt(p1_12cm*p1_12cm+M2sq), -p1mu.x(), -p1mu.y(), -p1mu.z()};
        // boost p1, p1 back to (1+2+3) CoM Frame
        p1mu = p1mu.boost_back(v12[0], v12[1], v12[2]); // p1
        p2mu = p2mu.boost_back(v12[0], v12[1], v12[2]); // p2
        // rotate p1, p2, p3(k) by phi
        double phi = Srandom::dist_phi(Srandom::gen);
        p1mu = p1mu.rotate_around_pz(phi);
        p2mu = p2mu.rotate_around_pz(phi);
        kmu = kmu.rotate_around_pz(phi);
        FS.clear();
        FS.push_back(p1mu);
        FS.push_back(p2mu);
        FS.push_back(kmu);
        return true;
    } 
    else{        
        LOG_INFO << "sample failed in 2->3 X";
	FS.clear();
        pids.clear();
        return false;
    }
}

/*****************************************************************/
/*******************find max of dX/dPS ***************************/
/*****************************************************************/
/*------------------Implementation for 2 -> 2--------------------*/
template<>
scalar Xsection<HS2PP, 2, double(*)(const double, void*)>::
    find_max(std::vector<double> parameters){
    double lnsqrts = parameters[0], temp = parameters[1];
    double sqrts = std::exp(lnsqrts);
    double s = std::pow(sqrts,2);
    double mA = _IS_masses[0], mB = _IS_masses[1];
    double m1 = _FS_masses[0], m2 = _FS_masses[1];
    double mD2 = t_channel_mD2->get_mD2(temp);
    double tcut = isPairProduction(_process_id)?
                -mD2/16.:-cut*mD2;
    double smin = std::max(std::pow(std::sqrt(mA*mA-tcut/2.)
                                   +std::sqrt(mB*mB-tcut/2.), 2),
                 std::pow(m1+m2,2) );
    s = std::max(smin*1.2, s);
    sqrts = std::sqrt(s);
    double pAcm = std::sqrt(
                  (std::pow(sqrts-mA,2)-mB*mB)
                * (std::pow(sqrts+mA,2)-mB*mB)
                  /4./s);
    double p1cm = std::sqrt(
                  (std::pow(sqrts-m1,2)-m2*m2)
                * (std::pow(sqrts+m1,2)-m2*m2)
                  /4./s);
    double tm2 = std::pow((mA*mA-m1*m1-mB*mB+m2*m2)/2./sqrts, 2);
    double tmin0 = tm2 - std::pow(pAcm+p1cm, 2);
    double tmax0 = tm2 - std::pow(pAcm-p1cm, 2);
    double tmin = tmin0-tcut;
    double tmax = tmax0+tcut;
    double cosmax = 1. - (tmax0-tmax)/(2.*pAcm*p1cm);
    double cosmin = 1. - (tmax0-tmin)/(2.*pAcm*p1cm);
    double etamax = .5*std::log((1.+cosmax)/(1.-cosmax));
    double etamin = .5*std::log((1.+cosmin)/(1.-cosmin));


    // Define a lamda function, which transforms the integral
    // varaible for more efficient integration
    auto minus_dXdeta = [s, temp, tmax, mA, mB, m1, m2, 
                        pAcm, p1cm, tmax0, this](double eta) {
        double params[7] = {s, tmax, temp, mA, mB, m1, m2};
        double costheta = std::tanh(eta);
        double t = tmax0 - 2.*pAcm*p1cm*(1.-costheta);
        double Jacobian = 2.*pAcm*p1cm*(1.-costheta*costheta);
        return -this->_f(t, params)*Jacobian;
    };
    double res = -minimize_1d(minus_dXdeta, {etamin, etamax}, 1e-8, 100, 1000);
    //LOG_INFO << sqrts << " " << temp << " " << res << " " << tmin << " " << tmax;
    return scalar{std::log(1.5*res)};
}

/*------------------Implementation for 2 -> 3--------------------*/
template<>
scalar Xsection<HS2PPP, 2, double(*)(const double*, void*)>::
    find_max(std::vector<double> parameters){
    double lnsqrts = parameters[0], temp = parameters[1];
    double mD2 = t_channel_mD2->get_mD2(temp);
    double tcut = is2to2_nlo(_process_id)?-mD2:-cut*mD2;
    double mA = _IS_masses[0], mB = _IS_masses[1];
    double m1 = _FS_masses[0], m2 = _FS_masses[1], m3 = _FS_masses[2];
    double sqrts = std::exp(lnsqrts);
    double s = std::pow(sqrts,2);
    double smin = std::max(
     std::pow(std::sqrt(mA*mA-tcut/4.)+std::sqrt(mB*mB-tcut/4.), 2), 
     std::pow(std::sqrt(m1*m1-tcut/9.)
             +std::sqrt(m2*m2-tcut/9.)+std::sqrt(m3*m3-tcut/9.), 2)
);
    s = std::max(2*smin, s);
    sqrts = std::sqrt(s);

    double pAcm = std::sqrt(
                  (std::pow(sqrts-mA,2)-mB*mB)
                * (std::pow(sqrts+mA,2)-mB*mB)
                  /4./s);
    //       pT3_sq 
    //       y3
    //       costheta12 // in 1+2 com frame
    //       phi12      // in 1+2 com frame
    // x0 = log(1+kt^2/mD2), 
    // x1 = y3/log(sqrts/kt)
    // x2 = log(1+2*(1-costheta12)*pAcm^2/mD2) ~ log(1+qt_cm12^2/mD2)
    // x3 = phi12
    // Jacobian = d(kT2, y3, c12, phi12)/d(x0, x1, x3, x3) = (kT2+mD2)*(c12+mD/2/pA/pA)*log(sqrts/kt)
    double mD2_2pA2 = mD2 / 2./pAcm/pAcm;
    double xmin[4] = {0., -1., 0., 0.};
    double xmax[4] = {std::log(1.+pAcm*pAcm/mD2), 1., std::log(1+2./mD2_2pA2), c2pi};
    auto dXdPS = [s, tcut, temp, mA, mB, m1, m2, m3, mD2, mD2_2pA2, pAcm, xmin, xmax, this](double * x_){
        if (x_[0]<xmin[0]||x_[0]>xmax[0]) return 0.;
        if (x_[1]<xmin[1]||x_[1]>xmax[1]) return 0.;
        if (x_[2]<xmin[2]||x_[2]>xmax[2]) return 0.;
        if (x_[3]<xmin[3]||x_[3]>xmax[3]) return 0.;
        double params[8] = {s, tcut, temp, mA, mB, m1, m2, m3};
        double kt2 = mD2*(std::exp(x_[0])-1);
        double ymax = std::acosh(pAcm/std::sqrt(kt2));
        double y3 = x_[1]*ymax, phi12 = x_[3]; 
        double costheta12 = 1.-mD2_2pA2*(std::exp(x_[2])-1);
        double xnew[4] = {kt2, y3, costheta12, phi12};
        double Jacobian = (1.-costheta12+ mD2_2pA2)*(kt2+mD2)*ymax;
        return this->_f(xnew, params)*Jacobian;
    };
    auto minus_dXdPS = [s, tcut, temp, mA, mB, m1, m2, m3, mD2, mD2_2pA2,pAcm, xmin, xmax, this](double * x_){
        if (x_[0]<xmin[0]||x_[0]>xmax[0]) return 0.;
        if (x_[1]<xmin[1]||x_[1]>xmax[1]) return 0.;
        if (x_[2]<xmin[2]||x_[2]>xmax[2]) return 0.;
        if (x_[3]<xmin[3]||x_[3]>xmax[3]) return 0.;
        double params[8] = {s, tcut, temp, mA, mB, m1, m2, m3};
        double kt2 = mD2*(std::exp(x_[0])-1);
        double ymax = std::acosh(pAcm/std::sqrt(kt2));
        double y3 = x_[1]*ymax, phi12 = x_[3]; 
        double costheta12 = 1.-mD2_2pA2*(std::exp(x_[2])-1);
        double xnew[4] = {kt2, y3, costheta12, phi12};
        double Jacobian = (1.-costheta12+ mD2_2pA2)*(kt2+mD2)*ymax;
        return -this->_f(xnew, params)*Jacobian;
    };
    // use MC_maximize to get into the vincinity ot the extrma
    auto startloc = MC_maximize(dXdPS, 4,
            {{xmin[0],xmax[0]}, 
             {xmin[1],xmax[1]},
             {xmin[2],xmax[2]}, 
             {xmin[3],xmax[3]}}, 1000);
    // use the best result of MC_maximize and determine the step of the simplex minimization method
    std::vector<double> step = {(xmax[0]-xmin[0])/20, (xmax[1]-xmin[1])/20., 
                                (xmax[2]-xmin[2])/20., (xmax[3]-xmin[3])/20.};
    for(int i=0; i<4; i++){
        double dx = std::min(xmax[i]-startloc[i], startloc[i]-xmin[i])/4.;
        step[i] = std::min(dx, step[i]);
    }
    // find the more precise maximum by the simplex method
    double val = -minimize_nd(minus_dXdPS, 4, startloc, step, 4000, 
   (xmax[0]-xmin[0])*(xmax[1]-xmin[1])*(xmax[2]-xmin[2])*(xmax[3]-xmin[3])
   /1e-12);
    return scalar{std::log(1.5*val)};

}


/*****************************************************************/
/*************************Integrate dX ***************************/
/*****************************************************************/
/*------------------Implementation for 2 -> 2--------------------*/
template<>
scalar Xsection<HS2PP, 2, double(*)(const double, void*)>::
        calculate_scalar(std::vector<double> parameters){
    double lnsqrts = parameters[0], temp = parameters[1];
    double sqrts = std::exp(lnsqrts);
    double s = std::pow(sqrts,2);
    double mA = _IS_masses[0], mB = _IS_masses[1];
    double m1 = _FS_masses[0], m2 = _FS_masses[1];
    if (s<=std::pow(m1+m2,2)||s<=std::pow(mA+mB,2)){
        return scalar{0.};
    }
    double pAcm = std::sqrt(
                  (std::pow(sqrts-mA,2)-mB*mB)
                * (std::pow(sqrts+mA,2)-mB*mB)
                  /4./s);
    double p1cm = std::sqrt(
                  (std::pow(sqrts-m1,2)-m2*m2)
                * (std::pow(sqrts+m1,2)-m2*m2)
                  /4./s);
    double mD2 = t_channel_mD2->get_mD2(temp);
    double tcut = isPairProduction(_process_id)?
                -mD2/16.:-cut*mD2;
    double tm2 = std::pow((mA*mA-m1*m1-mB*mB+m2*m2)/2./sqrts, 2);
    double tmin0 = tm2 - std::pow(pAcm+p1cm, 2);
    double tmax0 = tm2 - std::pow(pAcm-p1cm, 2);
    double tmin = tmin0-tcut;
    double tmax = tmax0+tcut;
    double cosmax = 1. - (tmax0-tmax)/(2.*pAcm*p1cm);
    double cosmin = 1. - (tmax0-tmin)/(2.*pAcm*p1cm);
    double etamax = .5*std::log((1.+cosmax)/(1.-cosmax));
    double etamin = .5*std::log((1.+cosmin)/(1.-cosmin));
    // Define a lamda function, which transforms the integral
    // varaible for more efficient integration
    auto dXdeta = [s, temp, tmax, mA, mB, m1, m2, 
                        pAcm, p1cm, tmax0, this](double eta) {
        double params[7] = {s, tmax, temp, mA, mB, m1, m2};
        double costheta = std::tanh(eta);
        double t = tmax0 - 2.*pAcm*p1cm*(1.-costheta);
        double Jacobian = 2.*pAcm*p1cm*(1.-costheta*costheta);
        return this->_f(t, params)*Jacobian;
    };
    if (tmin>=tmax) return scalar{0.};
    else {
        double error;
        double res = quad_1d(dXdeta, {etamin,etamax}, error);
        return scalar{res};
    }
}


/*------------------Implementation for 2 -> 3--------------------*/
template<>
scalar Xsection<HS2PPP, 2, double(*)(const double*, void*)>::
                calculate_scalar(std::vector<double> parameters){
    double lnsqrts = parameters[0], temp = parameters[1];
    double mD2 = t_channel_mD2->get_mD2(temp);
    double tcut = is2to2_nlo(_process_id)?-mD2:-cut*mD2;
    double mA = _IS_masses[0], mB = _IS_masses[1];
    double m1 = _FS_masses[0], m2 = _FS_masses[1], m3 = _FS_masses[2];
    double sqrts = std::exp(lnsqrts);
    double s = std::pow(sqrts,2);
    double smin = std::max(
     std::pow(std::sqrt(mA*mA-tcut/4.)+std::sqrt(mB*mB-tcut/4.), 2), 
     std::pow(std::sqrt(m1*m1-tcut/9.)
             +std::sqrt(m2*m2-tcut/9.)+std::sqrt(m3*m3-tcut/9.), 2)
);
    if (s<smin) return scalar{0.};
    sqrts = std::sqrt(s);
    double pAcm = std::sqrt(
                  (std::pow(sqrts-mA,2)-mB*mB)
                * (std::pow(sqrts+mA,2)-mB*mB)
                  /4./s);
    //       pT3_sq 
    //       y3
    //       costheta12 // in 1+2 com frame
    //       phi12      // in 1+2 com frame
    // x0 = log(1+kt^2/mD2), 
    // x1 = y3/log(sqrts/kt)
    // x2 = log(1+2*(1-costheta12)*pAcm^2/mD2) ~ log(1+qt_cm12^2/mD2)
    // x3 = phi12
    // Jacobian = d(kT2, y3, c12, phi12)/d(x0, x1, x3, x3) = (kT2+mD2)*(1-c12+mD/2/pA/pA)*log(sqrts/kt)
    double mD2_2pA2 = mD2 / (2.*pAcm*pAcm);
    auto dXdPS = [s, tcut, temp, mA, mB, m1, m2, m3, mD2, mD2_2pA2, pAcm, this](double * x_){
        double params[8] = {s, tcut, temp, mA, mB, m1, m2, m3};
        double kt2 = mD2*(std::exp(x_[0])-1);
        double ymax = std::acosh(pAcm/std::sqrt(kt2));
        double y3 = x_[1]*ymax, phi12 = x_[3]; 
        double costheta12 = 1.-mD2_2pA2*(std::exp(x_[2])-1);
        double xnew[4] = {kt2, y3, costheta12, phi12};
        double Jacobian = (1.-costheta12+ mD2_2pA2)*(kt2+mD2)*ymax;
        return this->_f(xnew, params)*Jacobian;
    };
    double xmin[4] = {0., -1.,  0., 0.};
    double xmax[4] = {std::log(1.+pAcm*pAcm/mD2), 1., std::log(1.+2./mD2_2pA2), c2pi};
    double error;
    double res = vegas(dXdPS, 4, xmin, xmax, error);
    return scalar{res};
}

/*****************************************************************/
/**************Integrate dX \Delta p^mu***************************/
/*****************************************************************/
/*------------------Default Implementation-----------------------*/
template<const char * str, size_t N, typename F>
fourvec Xsection<str, N, F>::calculate_fourvec(std::vector<double> parameters){
    return fourvec::unity();
}
/*------------------Implementation for 2 -> 2--------------------*/
template<>
fourvec Xsection<HS2PP, 2, double(*)(const double, void*)>::
        calculate_fourvec(std::vector<double> parameters){
    double lnsqrts = parameters[0], temp = parameters[1];
    double sqrts = std::exp(lnsqrts);
    double s = std::pow(sqrts,2);
    double mA = _IS_masses[0], mB = _IS_masses[1];
    double m1 = _FS_masses[0], m2 = _FS_masses[1];
    double pAcm = std::sqrt(
                  (std::pow(sqrts-mA,2)-mB*mB)
                * (std::pow(sqrts+mA,2)-mB*mB)
                  /4./s);
    double p1cm = std::sqrt(
                  (std::pow(sqrts-m1,2)-m2*m2)
                * (std::pow(sqrts+m1,2)-m2*m2)
                  /4./s);
    double tcut = -cut*t_channel_mD2->get_mD2(temp);
    double tm2 = std::pow((mA*mA-m1*m1-mB*mB+m2*m2)/2./sqrts, 2);
    double tmin0 = tm2 - std::pow(pAcm+p1cm, 2);
    double tmax0 = tm2 - std::pow(pAcm-p1cm, 2);
    double tmin = tmin0-tcut;
    double tmax = tmax0+tcut;
    double cosmax = 1. - (tmax0-tmax)/(2.*pAcm*p1cm);
    double cosmin = 1. - (tmax0-tmin)/(2.*pAcm*p1cm);
    double etamax = .5*std::log((1.+cosmax)/(1.-cosmax));
    double etamin = .5*std::log((1.+cosmin)/(1.-cosmin));
    // Define a lamda function, which transforms the integral
    // varaible for more efficient integration
    auto dXdeta_Pmu = [s, temp, tmax, mA, mB, m1, m2, 
                        pAcm, p1cm, tmax0, this](double eta) {
        double params[7] = {s, tmax, temp, mA, mB, m1, m2};
        double costheta = std::tanh(eta);
        double t = tmax0 - 2.*pAcm*p1cm*(1.-costheta);
        double Jacobian = 2.*pAcm*p1cm*(1.-costheta*costheta);
        double dpz = p1cm*costheta-pAcm;
        return this->_f(t, params)*dpz*Jacobian;
    };

    if (tmin>=tmax)  return fourvec{0., 0., 0., 0.};
    else {
        double error;
        double dpz = quad_1d(dXdeta_Pmu, {etamin, etamax}, error);
        return fourvec{0., 0., 0., dpz};
    }
}

/*****************************************************************/
/**************Integrate dX \Delta p^mu*p^nu *********************/
/*****************************************************************/
/*------------------Default Implementation-----------------------*/
template<const char * str, size_t N, typename F>
tensor Xsection<str, N, F>::calculate_tensor(std::vector<double> parameters){
    return tensor::unity();
}
/*------------------Implementation for 2 -> 2--------------------*/
template<>
tensor Xsection<HS2PP, 2, double(*)(const double, void*)>::
    calculate_tensor(std::vector<double> parameters){
        double lnsqrts = parameters[0], temp = parameters[1];
    double sqrts = std::exp(lnsqrts);
    double s = std::pow(sqrts,2);
    double mA = _IS_masses[0], mB = _IS_masses[1];
    double m1 = _FS_masses[0], m2 = _FS_masses[1];
    double pAcm = std::sqrt(
                  (std::pow(sqrts-mA,2)-mB*mB)
                * (std::pow(sqrts+mA,2)-mB*mB)
                  /4./s);
    double p1cm = std::sqrt(
                  (std::pow(sqrts-m1,2)-m2*m2)
                * (std::pow(sqrts+m1,2)-m2*m2)
                  /4./s);
    double tcut = -cut*t_channel_mD2->get_mD2(temp);
    double tm2 = std::pow((mA*mA-m1*m1-mB*mB+m2*m2)/2./sqrts, 2);
    double tmin0 = tm2 - std::pow(pAcm+p1cm, 2);
    double tmax0 = tm2 - std::pow(pAcm-p1cm, 2);
    double tmin = tmin0-tcut;
    double tmax = tmax0+tcut;
    double cosmax = 1. - (tmax0-tmax)/(2.*pAcm*p1cm);
    double cosmin = 1. - (tmax0-tmin)/(2.*pAcm*p1cm);
    double etamax = .5*std::log((1.+cosmax)/(1.-cosmax));
    double etamin = .5*std::log((1.+cosmin)/(1.-cosmin));
    // Define a lamda function, which transforms the integral
    // varaible for more efficient integration
    auto dXdeta_PzPz = [s, temp, tmax, mA, mB, m1, m2, 
                        pAcm, p1cm, tmax0, this](double eta) {
        double params[7] = {s, tmax, temp, mA, mB, m1, m2};
        double costheta = std::tanh(eta);
        double t = tmax0 - 2.*pAcm*p1cm*(1.-costheta);
        double Jacobian = 2.*pAcm*p1cm*(1.-costheta*costheta);
        double dpz = p1cm*costheta-pAcm;
        return this->_f(t, params)*dpz*dpz*Jacobian;
    };

    auto dXdeta_PTPT = [s, temp, tmax, mA, mB, m1, m2, 
                        pAcm, p1cm, tmax0, this](double eta) {
        double params[7] = {s, tmax, temp, mA, mB, m1, m2};
        double costheta = std::tanh(eta);
        double t = tmax0 - 2.*pAcm*p1cm*(1.-costheta);
        double Jacobian = 2.*pAcm*p1cm*(1.-costheta*costheta);
        double dpt2 = p1cm*p1cm*(1.-costheta*costheta);
        return this->_f(t, params)*dpt2*Jacobian;
    };


    double dpzdpz, dptdpt;
    if (tmin>tmax) {
        dpzdpz=0.;
        dptdpt=0.;
    }
    else{
        double error;
        dpzdpz = quad_1d(dXdeta_PzPz, {etamin, etamax}, error);
        dptdpt = quad_1d(dXdeta_PTPT, {etamin, etamax}, error);
    }
    return tensor{0.,     0.,         0.,         0.,
                  0.,     dptdpt/2.,     0.,         0.,
                  0.,     0.,         dptdpt/2.,    0.,
                  0.,     0.,         0.,         dpzdpz};
}
// instance:
template class Xsection<HS2PP, 2, double(*)(const double, void*)>;
template class Xsection<HS2PPP, 2, double(*)(const double*, void*)>;
