#include "Rate.h"
#include "integrator.h"
#include "sampler.h"
#include "minimizer.h"
#include "random.h"
#include "approx_functions.h"
#include "matrix_elements.h"
#include "Langevin.h"
#include <iostream>

template <>
Rate<HS2PP, 2, 2, double(*)(const double, void *)>::
    Rate(std::string Name, std::string configfile, double(*f)(const double, void *)):
StochasticBase<2>(Name+"/rate", configfile),
X(std::make_shared<Xsection<HS2PP, 2, double(*)(const double, void *)>>(Name, configfile, f) )
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

    _degen = tree.get<double>("degeneracy");
    _active = (tree.get<std::string>("<xmlattr>.status")=="active")?true:false;

}

template <>
Rate<HS2PPP, 2, 2, double(*)(const double*, void *)>::
    Rate(std::string Name, std::string configfile, double(*f)(const double*, void *)):
StochasticBase<2>(Name+"/rate", configfile),
X(std::make_shared<Xsection<HS2PPP, 2, double(*)(const double*, void *)>>(Name, configfile, f) )
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
    _degen = tree.get<double>("degeneracy");
    _active = (tree.get<std::string>("<xmlattr>.status")=="active")?true:false;
}


/*****************************************************************/
/*************************Sample dR ******************************/
/*****************************************************************/
/*------------------Implementation for 2 -> 2--------------------*/
template <>
bool Rate<HS2PP, 2, 2, double(*)(const double, void *)>::
        sample(std::vector<double> parameters,
                        int incoming_hard_pid,
			std::vector<fourvec> & final_states,
                        std::vector<int> & pids ){
    double lnEA = parameters[0], T = parameters[1];
    double EA = std::exp(lnEA);
    double mA = _IS_masses[0], mB = _IS_masses[1];
    double m1 = _FS_masses[0], m2 = _FS_masses[1];
    double pA = std::sqrt(EA*EA-mA*mA);
    double mD2 = t_channel_mD2->get_mD2(T);
    double tcut = isPairProduction(_process_id)?
                -mD2/16.:-cut*mD2;
    double smin = std::max(
     std::pow(std::sqrt(mA*mA-tcut/2.)+std::sqrt(mB*mB-tcut/2.), 2), 
     std::pow(std::sqrt(m1*m1-tcut/2.)+std::sqrt(m2*m2-tcut/2.), 2) );
    auto dR_dxdy = [EA, T, smin, pA, mA, mB, this](const double * x){
        double pB = x[0], costheta = x[1];
        double EB = std::sqrt(pB*pB+mB*mB);
        double s = mA*mA + mB*mB + 2.*(EA*EB-pA*pB*costheta);
        if (s<smin*1.5||std::abs(costheta)>=1.) return 0.;
        double flux = 2.*std::sqrt((s-std::pow(mA-mB,2))
                                  *(s-std::pow(mA+mB,2)));
        double Xtot = this->X->GetZeroM({.5*std::log(s),T}).s;
        return pB*pB/EB/EA*std::exp(-EB/T)*flux*Xtot/c16pi2;
    };
    double pBmin = (smin-mA*mA-mB*mB)/(2.*(EA*EA-mA*mA)) 
           * (std::sqrt(1.-4*mA*mA*mB*mB/std::pow(smin-mA*mA-mB*mB,2))*EA-pA);
    double pBmax = pBmin+8.*T;
    bool status = true;
    double fmax = StochasticBase<2>::GetFmax(parameters).s;
    auto res = sample_nd(dR_dxdy, 2, 
                      {{pBmin, pBmax}, {-1., 1.}},
                      fmax, status);
    if (status == true){
        double pB = res[0];
        double EB = std::sqrt(pB*pB+mB*mB);
        double costheta = res[1];
        double s = mA*mA + mB*mB + 2.*(EA*EB-pA*pB*costheta);
        double lnsqrts = .5*std::log(s);
        double sintheta = std::sqrt(1. - costheta*costheta);
        bool statusX = X->sample({lnsqrts, T}, 
                                  incoming_hard_pid, final_states, pids);
        if (!statusX) {
            LOG_INFO << "sample rate failed--1, 2->2";
            final_states.clear();
            pids.clear();
            return statusX;
        }
        // give incoming partilce a random phi angle
        double phi = Srandom::dist_phi(Srandom::gen);
        // com velocity
        double vcom[3] = { EB*sintheta*std::cos(phi)/(EB+EA),
                           EB*sintheta*std::sin(phi)/(EB+EA),
                          (EB*costheta+pA)/(EB+EA)    };
        /*
        FS now is in Z-oriented CoM frame
        1) FS.rotate_back
        2) FS.boost_back
        */
        fourvec p1{EA, 0, 0, pA};
        auto p1com = p1.boost_to(vcom[0], vcom[1], vcom[2]);
        for(auto & p: final_states){
            p = p.rotate_back(p1com);
            p = p.boost_back(vcom[0], vcom[1], vcom[2]);
        }
        return statusX;
    }
    else{
        LOG_INFO << "sample rate failed--2, 2->2";
        LOG_INFO <<"prcid = "<< _process_id<< ", EA = " << EA << ", T = " << T;
        final_states.clear();
        pids.clear();
        return false;
    }
}




/*------------------Implementation for 2 -> 3--------------------*/
template <>
bool Rate<HS2PPP, 2, 2, double(*)(const double*, void *)>::
        sample(std::vector<double> parameters,
                        int incoming_hard_pid,
			std::vector<fourvec> & final_states,
                        std::vector<int> & pids ){
    double lnEA = parameters[0], T = parameters[1];
    double EA = std::exp(lnEA);
    double mA = _IS_masses[0], mB = _IS_masses[1];
    double m1 = _FS_masses[0], m2 = _FS_masses[1], m3 = _FS_masses[2];
    double pA = std::sqrt(EA*EA-mA*mA);
    double tcut = -cut*t_channel_mD2->get_mD2(T);
    double smin = std::max(
     std::pow(std::sqrt(mA*mA-tcut/4.)+std::sqrt(mB*mB-tcut/4.), 2), 
     std::pow(std::sqrt(m1*m1-tcut/9.)+std::sqrt(m2*m2-tcut/9.)
            +std::sqrt(m3*m3-tcut/9.), 2) );

    auto dR_dxdy = [EA, T, smin, pA, mA, mB, this](const double * x){
        double pB = x[0], costheta = x[1];
        double EB = std::sqrt(pB*pB+mB*mB);
        double s = mA*mA + mB*mB + 2.*(EA*EB-pA*pB*costheta);
        if (s<smin||std::abs(costheta)>=1.) return 0.;
        double flux = 2.*std::sqrt((s-std::pow(mA-mB,2))
                                  *(s-std::pow(mA+mB,2)));
        double Xtot = this->X->GetZeroM({.5*std::log(s),T}).s;
        return pB*pB/EB/EA*std::exp(-EB/T)*flux*Xtot/c16pi2;
    };
    double pBmin = (smin-mA*mA-mB*mB)/(2.*(EA*EA-mA*mA)) 
           * (std::sqrt(1.-4*mA*mA*mB*mB/std::pow(smin-mA*mA-mB*mB,2))*EA-pA);
    double pBmax = pBmin+8.*T;
    bool status = true;
    double fmax = StochasticBase<2>::GetFmax(parameters).s;
    auto res = sample_nd(dR_dxdy, 2, 
                      {{pBmin, pBmax}, {-1., 1.}},
                      fmax, status);
    if (status==true){
        double pB = res[0];
        double EB = std::sqrt(pB*pB+mB*mB);
        double costheta = res[1];
        double s = mA*mA + mB*mB + 2.*(EA*EB-pA*pB*costheta);
        double lnsqrts = .5*std::log(s);
        double sintheta = std::sqrt(1. - costheta*costheta);
        bool statusX = X->sample({lnsqrts, T}, 
                                  incoming_hard_pid, final_states, pids);
        if (!statusX) {
            LOG_INFO << "sample rate failed--1, 2->3";
            final_states.clear();
            pids.clear();
            return statusX;
        }
        // give incoming partilce a random phi angle
        double phi = Srandom::dist_phi(Srandom::gen);
        // com velocity
        double vcom[3] = { EB*sintheta*std::cos(phi)/(EB+EA),
                           EB*sintheta*std::sin(phi)/(EB+EA),
                          (EB*costheta+pA)/(EB+EA)    };
        /*
        FS now is in Z-oriented CoM frame
        1) FS.rotate_back
        2) FS.boost_back
        */
        fourvec p1{EA, 0, 0, pA};
        auto p1com = p1.boost_to(vcom[0], vcom[1], vcom[2]);
        for(auto & p: final_states){
            p = p.rotate_back(p1com);
            p = p.boost_back(vcom[0], vcom[1], vcom[2]);
        }
        return statusX;
    }
    else{
        LOG_INFO << "sample rate failed 2->3";
        final_states.clear();
        pids.clear();
        return false;
    }
}

/*****************************************************************/
/*************************Find dR_max ****************************/
/*****************************************************************/
/*------------------Implementation for 2 -> 2--------------------*/
template <>
scalar Rate<HS2PP, 2, 2, double(*)(const double, void*)>::
        find_max(std::vector<double> parameters){
    double lnEA = parameters[0], T = parameters[1];
    double EA = std::exp(lnEA);
    double mA = _IS_masses[0], mB = _IS_masses[1];
    double m1 = _FS_masses[0], m2 = _FS_masses[1];
    double pA = std::sqrt(EA*EA-mA*mA);
    double mD2 = t_channel_mD2->get_mD2(T);
    double tcut = isPairProduction(_process_id)?
                -mD2/16.:-cut*mD2;
    double smin = std::max(
     std::pow(std::sqrt(mA*mA-tcut/2.)+std::sqrt(mB*mB-tcut/2.), 2), 
     std::pow(std::sqrt(m1*m1-tcut/2.)+std::sqrt(m2*m2-tcut/2.), 2) );
    auto minus_dR_dxdy = [EA, T, smin, pA, mA, mB, this](const double * x){
        double pB = x[0], costheta = x[1];
        double EB = std::sqrt(pB*pB+mB*mB);
        double s = mA*mA + mB*mB + 2.*(EA*EB-pA*pB*costheta);
        if (s<smin*1.5||std::abs(costheta)>=1.) return 0.;
        double flux = 2.*std::sqrt((s-std::pow(mA-mB,2))
                                  *(s-std::pow(mA+mB,2)));
        double Xtot = this->X->GetZeroM({.5*std::log(s),T}).s;
        return -pB*pB/EB/EA*std::exp(-EB/T)*flux*Xtot/c16pi2;
    };
    double pBmin = (smin-mA*mA-mB*mB)/(2.*(EA*EA-mA*mA)) 
           * (std::sqrt(1.-4*mA*mA*mB*mB/std::pow(smin-mA*mA-mB*mB,2))*EA-pA);
    double pBmax = pBmin+8.*T;
    auto val = -minimize_nd(minus_dR_dxdy, 2, 
                        {(pBmin+pBmax)/2., -1.}, {(pBmax-pBmin)/20., 0.1}, 
                          1000, 1e-12);
    return scalar{1.5*val};
}


/*------------------Implementation for 2 -> 3--------------------*/
template <>
scalar Rate<HS2PPP, 2, 2, double(*)(const double*, void*)>::
        find_max(std::vector<double> parameters){
    double lnEA = parameters[0], T = parameters[1];
    double EA = std::exp(lnEA);
    double mA = _IS_masses[0], mB = _IS_masses[1];
    double m1 = _FS_masses[0], m2 = _FS_masses[1], m3 = _FS_masses[2];
    double pA = std::sqrt(EA*EA-mA*mA);
    double tcut = -cut*t_channel_mD2->get_mD2(T);
    double smin = std::max(
     std::pow(std::sqrt(mA*mA-tcut/4.)+std::sqrt(mB*mB-tcut/4.), 2), 
     std::pow(std::sqrt(m1*m1-tcut/9.)+std::sqrt(m2*m2-tcut/9.)
             +std::sqrt(m3*m3-tcut/9.), 2) );
    auto minus_dR_dxdy = [EA, T, smin, pA, mA, mB, this](const double * x){
        double pB = x[0], costheta = x[1];
        double EB = std::sqrt(pB*pB+mB*mB);
        double s = mA*mA + mB*mB + 2.*(EA*EB-pA*pB*costheta);
        if (s<smin||std::abs(costheta)>=1.) return 0.;
        double flux = 2.*std::sqrt((s-std::pow(mA-mB,2))
                                  *(s-std::pow(mA+mB,2)));
        double Xtot = this->X->GetZeroM({.5*std::log(s),T}).s;
        return -pB*pB/EB/EA*std::exp(-EB/T)*flux*Xtot/c16pi2;
    };
    double pBmin = (smin-mA*mA-mB*mB)/(2.*(EA*EA-mA*mA)) 
           * (std::sqrt(1.-4*mA*mA*mB*mB/std::pow(smin-mA*mA-mB*mB,2))*EA-pA);
    double pBmax = pBmin+8.*T;
    auto val = -minimize_nd(minus_dR_dxdy, 2, 
                        {(pBmin+pBmax)/2., -1.}, {(pBmax-pBmin)/10., 0.2}, 
                          1000, 1e-12);
    return scalar{1.5*val};
}

/*****************************************************************/
/*************************Integrate dR ***************************/
/*****************************************************************/
/*------------------Implementation for 2 -> 2--------------------*/
template <>
scalar Rate<HS2PP, 2, 2, double(*)(const double, void*)>::
        calculate_scalar(std::vector<double> parameters){
    double lnEA = parameters[0], T = parameters[1];
    double EA = std::exp(lnEA);
    double mA = _IS_masses[0], mB = _IS_masses[1];
    double m1 = _FS_masses[0], m2 = _FS_masses[1];
    double pA = std::sqrt(EA*EA-mA*mA);
    double tcut = -cut*t_channel_mD2->get_mD2(T);
    double smin = std::max(
     std::pow(std::sqrt(mA*mA-tcut/4.)+std::sqrt(mB*mB-tcut/4.), 2), 
     std::pow(std::sqrt(m1*m1-tcut/4.)+std::sqrt(m2*m2-tcut/4.), 2) );
    auto dR_dxdy = [EA, T, smin, pA, mA, mB, this](const double * x){
        double pB = x[0], costheta = x[1];
        double EB = std::sqrt(pB*pB+mB*mB);
        double s = mA*mA + mB*mB + 2.*(EA*EB-pA*pB*costheta);
        double result = 0.;
        if (s<smin||std::abs(costheta)>=1.) result = 0.;
        else {
            double flux = 2.*std::sqrt((s-std::pow(mA-mB,2))
                                      *(s-std::pow(mA+mB,2)));
            double Xtot = this->X->GetZeroM({.5*std::log(s),T}).s;
            result = pB*pB/EB/EA*std::exp(-EB/T)*flux*Xtot/c16pi2;
        }
        std::vector<double> res{result};
        return res;
    };
    double pBmin = (smin-mA*mA-mB*mB)/(2.*(EA*EA-mA*mA)) 
           * (std::sqrt(1.-4*mA*mA*mB*mB/std::pow(smin-mA*mA-mB*mB,2))*EA-pA);
    double pBmax = pBmin+8.*T;
    double xmin[2] = {pBmin, -1.};
    double xmax[2] = {pBmax,1.};
    double err;
    auto val = quad_nd(dR_dxdy, 2, 1, xmin, xmax, err)[0] * _degen;
    return scalar{val};
}

/*------------------Implementation for 2 -> 3--------------------*/
template <>
scalar Rate<HS2PPP, 2, 2, double(*)(const double*, void*)>::
        calculate_scalar(std::vector<double> parameters){
    double lnEA = parameters[0], T = parameters[1];
    double EA = std::exp(lnEA);
    double mA = _IS_masses[0], mB = _IS_masses[1];
    double m1 = _FS_masses[0], m2 = _FS_masses[1], m3 = _FS_masses[2];
    double pA = std::sqrt(EA*EA-mA*mA);
    double tcut = -cut*t_channel_mD2->get_mD2(T);
    double smin = std::max(
     std::pow(std::sqrt(mA*mA-tcut/4.)+std::sqrt(mB*mB-tcut/4.), 2), 
     std::pow(std::sqrt(m1*m1-tcut/9.)+std::sqrt(m2*m2-tcut/9.)
             +std::sqrt(m3*m3-tcut/9.), 2) );
    auto dR_dxdy = [EA, T, smin, pA, mA, mB, this](const double * x){
        double pB = x[0], costheta = x[1];
        double EB = std::sqrt(pB*pB+mB*mB);
        double s = mA*mA + mB*mB + 2.*(EA*EB-pA*pB*costheta);
        double result = 0.;
        if (s<smin||std::abs(costheta)>=1.) result = 0.;
        else {
            double flux = 2.*std::sqrt((s-std::pow(mA-mB,2))
                                      *(s-std::pow(mA+mB,2)));
            double Xtot = this->X->GetZeroM({.5*std::log(s),T}).s;
            result = pB*pB/EB/EA*std::exp(-EB/T)*flux*Xtot/c16pi2;
        }
        std::vector<double> res{result};
        return res;
    };
    double pBmin = (smin-mA*mA-mB*mB)/(2.*(EA*EA-mA*mA)) 
           * (std::sqrt(1.-4*mA*mA*mB*mB/std::pow(smin-mA*mA-mB*mB,2))*EA-pA);
    double pBmax = pBmin+8.*T;
    double xmin[2] = {pBmin, -1.};
    double xmax[2] = {pBmax,1.};
    double err;
    auto val = quad_nd(dR_dxdy, 2, 1, xmin, xmax, err)[0] * _degen;
    return scalar{val};
}

/*****************************************************************/
/*******************Integrate Delta_p^mu*dR **********************/
/*****************************************************************/
/*------------------No need to do this generally-----------------*/
template <char const *str, size_t N1, size_t N2, typename F>
fourvec Rate<str, N1, N2, F>::calculate_fourvec(std::vector<double> parameters){
    return fourvec::unity();
}
/*------------------Implementation for 2 -> 2--------------------*/
template <>
fourvec Rate<HS2PP, 2, 2, double(*)(const double, void*)>::
        calculate_fourvec(std::vector<double> parameters){
    double lnEA = parameters[0], T = parameters[1];
    double EA = std::exp(lnEA);
    double mA = _IS_masses[0], mB = _IS_masses[1];
    double m1 = _FS_masses[0], m2 = _FS_masses[1];
    double pA = std::sqrt(EA*EA-mA*mA);
    double tcut = -cut*t_channel_mD2->get_mD2(T);
    double smin = std::max(
     std::pow(std::sqrt(mA*mA-tcut/4.)+std::sqrt(mB*mB-tcut/4.), 2), 
     std::pow(std::sqrt(m1*m1-tcut/4.)+std::sqrt(m2*m2-tcut/4.), 2) );
    auto code = [EA, T, smin, pA, mA, mB, this](const double * x){
        double pB = x[0], costheta = x[1];
        double EB = std::sqrt(pB*pB+mB*mB);
        double s = mA*mA + mB*mB + 2.*(EA*EB-pA*pB*costheta);
        if (s<smin||std::abs(costheta)>=1.) {
            std::vector<double> res{0., 0.};
            return res;
        }
        else {
            double flux = 2.*std::sqrt((s-std::pow(mA-mB,2))
                                      *(s-std::pow(mA+mB,2)));
            double lnsqrts = .5*std::log(s);
            // A vector in p1z(com)-oriented com frame
            auto fmu0 = this->X->GetFirstM({lnsqrts,T});
            // rotate it back from p1z(com)-oriented com frame
            fourvec p1{EA, 0, 0, pA};
            double sintheta = std::sqrt(1. - costheta*costheta);
            double vcom[3] = {EB*sintheta/(EA+EB), 0., 
                             (EB*costheta+pA)/(EA+EB)};
            auto fmu1 = fmu0.rotate_back(
                   p1.boost_to(vcom[0], vcom[1], vcom[2]));
            // boost back to the matter frame
            auto fmu2 = fmu1.boost_back(vcom[0], vcom[1], vcom[2]);
            double C = pB*pB/EB/EA*std::exp(-EB/T)*flux/c16pi2;
            std::vector<double> res{C*fmu2.t(), C*fmu2.z()};
            return res;
        }
    };
    double pBmin = (smin-mA*mA-mB*mB)/(2.*(EA*EA-mA*mA)) 
           * (std::sqrt(1.-4*mA*mA*mB*mB/std::pow(smin-mA*mA-mB*mB,2))*EA-pA);
    double pBmax = pBmin+8.*T;
    double xmin[2] = {pBmin, -1.};
    double xmax[2] = {pBmax, 1.};
    double err;
    auto val = quad_nd(code, 2, 2, xmin, xmax, err);
    return fourvec{_degen*val[0], 0.0, 0.0, _degen*val[1]};
}

/*****************************************************************/
/***********Integrate Delta_p^mu*Delta_p^mu*dR *******************/
/*****************************************************************/
/*------------------No need to do this generally-----------------*/
template <char const *str, size_t N1, size_t N2, typename F>
tensor Rate<str, N1, N2, F>::calculate_tensor(std::vector<double> parameters){
    return tensor::unity();
}
/*------------------Implementation for 2 -> 2--------------------*/
template <>
tensor Rate<HS2PP, 2, 2, double(*)(const double, void*)>::
        calculate_tensor(std::vector<double> parameters){
    double lnEA = parameters[0], T = parameters[1];
    double EA = std::exp(lnEA);
    double mA = _IS_masses[0], mB = _IS_masses[1];
    double m1 = _FS_masses[0], m2 = _FS_masses[1];
    double pA = std::sqrt(EA*EA-mA*mA);
    double tcut = -cut*t_channel_mD2->get_mD2(T);
    double smin = std::max(
     std::pow(std::sqrt(mA*mA-tcut/4.)+std::sqrt(mB*mB-tcut/4.), 2), 
     std::pow(std::sqrt(m1*m1-tcut/4.)+std::sqrt(m2*m2-tcut/4.), 2) );
    auto code = [EA, T, smin, pA, mA, mB, this](const double * x){
        double pB = x[0], costheta = x[1], phi = x[2];
        double EB = std::sqrt(pB*pB+mB*mB);
        double s = mA*mA + mB*mB + 2.*(EA*EB-pA*pB*costheta);
        if (s<smin||costheta>=1.||costheta<=-1.||phi<=0.||phi>c2pi) {
            std::vector<double> res{0., 0., 0., 0.};
            return res;
        }
        else{
            double flux = 2.*std::sqrt((s-std::pow(mA-mB,2))
                                      *(s-std::pow(mA+mB,2)));
            double lnsqrts = .5*std::log(s);
            double sintheta = std::sqrt(1. - costheta*costheta);
            double vcom[3] = {EB*sintheta/(EB+EA)*std::cos(phi), 
                    EB*sintheta/(EB+EA)*std::sin(phi),
                   (EB*costheta+pA)/(EB+EA)};
            fourvec p1{EA, 0, 0, pA};
            // A vector in p1z(com)-oriented com frame
            auto fmunu0 = this->X->GetSecondM({lnsqrts,T});
            // rotate it back from p1z(com)-oriented com frame
            auto fmunu1 = fmunu0.rotate_back(
                           p1.boost_to(vcom[0], vcom[1], vcom[2]));
            // boost back to the matter frame
            auto fmunu2 = fmunu1.boost_back(vcom[0], vcom[1], vcom[2]);
            double C = pB*pB/EB/EA*std::exp(-EB/T)*flux/c32pi3;
            // Set tranverse to zero due to azimuthal symmetry;
            std::vector<double> res{C*fmunu2.T[0][0], C*fmunu2.T[1][1],
                                C*fmunu2.T[2][2], C*fmunu2.T[3][3]};
            return res;
        }
    };
    double pBmin = (smin-mA*mA-mB*mB)/(2.*(EA*EA-mA*mA)) 
           * (std::sqrt(1.-4*mA*mA*mB*mB/std::pow(smin-mA*mA-mB*mB,2))*EA-pA);
    double pBmax = pBmin+8.*T;
    double xmin[3] = {pBmin, -1., 0.};
    double xmax[3] = {pBmax, 1., c2pi};
    double err;
    auto val = quad_nd(code, 3, 4, xmin, xmax, err);
    return tensor{_degen*val[0], 0., 0., 0.,
                  0., _degen*val[1], 0., 0.,
                  0., 0., _degen*val[2], 0.,
                  0., 0., 0., _degen*val[3]};
}


//Diff-1-to-2 Rate
template<>
EffRate12<2, double(*)(const double*, void *)>::
    EffRate12(std::string Name, std::string configfile, double(*f)(const double *, void *)):
StochasticBase<2>(Name+"/rate", configfile),
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
   auto tree = config.get_child(model_name + "." + process_name);
   _degen = tree.get<double>("degeneracy");
   _process_id = get_process_info(process_name, _IS_masses, _FS_masses, 
                                   _IS_types, _FS_types);
   _active = (tree.get<std::string>("<xmlattr>.status")=="active")?true:false;
}

// Sample Final states
template<>
bool EffRate12<2, double(*)(const double*, void *)>::
    sample(std::vector<double> parameters,
                        int incoming_hard_pid,
			std::vector<fourvec> & final_states,
                        std::vector<int> & pids ){
    double lnEA = parameters[0], T = parameters[1];
    double EA = std::exp(lnEA);
    double mA = _IS_masses[0];
    double m1 = _FS_masses[0], m2 = _FS_masses[1];
    EA = std::max(std::max(EA, mA), m1+m2);
    double pA = std::sqrt(EA*EA - mA*mA);
    double lnxmin = std::log((.01+m1)/(EA+pA));
    double lnxmax = std::log(1.-m2/(EA+pA));
    double lnsinmin = std::min(std::log(.01/EA), -9.);
    double lnsinmax = 0.;
    auto dR_dlnxdlny = [EA, T, mA, m1, m2, lnxmin, lnxmax, lnsinmin, lnsinmax, this](const double *x){
        if (x[0]<lnxmin||x[0]>lnxmax||x[1]<lnsinmin||x[1]>lnsinmax) return 0.;
        double X[2] = {std::exp(x[0]), std::exp(x[1])};
        double Jacobian = X[0]*X[1];
        double params[5] = {EA, T, mA, m1, m2};
        return this->_f(X, params)*Jacobian;
    };
    bool status=true;
    double fmax = std::exp(StochasticBase<2>::GetFmax(parameters).s);
    auto res = sample_nd(dR_dlnxdlny, 2, 
                        {{lnxmin ,lnxmax}, {lnsinmin, lnsinmax}}, 
                         fmax, status);
    if (status==true){
        double x = std::exp(res[0]); // k+/p+
        double sintheta = std::exp(res[1]); // sintheta = kT/kabs > 0
        double tantheta = sintheta/std::sqrt(1.-sintheta*sintheta);
        double kplus = x*(EA+pA);
        double kT = std::sqrt((kplus/sintheta+m1)*(kplus/sintheta-m1))
                  - kplus/tantheta;
        double kabs = kT/sintheta;
        double phi = Srandom::dist_phi(Srandom::gen);
        double kx = kT*std::cos(phi), ky = kT*std::sin(phi);
        double kz = kT/tantheta;
        double k0 = (kplus - kz);
        double qz = pA-kz;
        double q0 = std::sqrt(qz*qz + kT*kT + m1*m1);
        assign_1to2_pid(_process_id, incoming_hard_pid, pids);
        final_states.resize(2);
        final_states[0] = fourvec{q0, -kx, -ky, qz};
        final_states[1] = fourvec{k0, kx, ky, kz};
        return true;
    }
    else{
        LOG_INFO << "sample rate failed, 1->2";
        final_states.clear();
        pids.clear();
        return false;
    }
}

// Find the max of dR/dx/dy
template <>
scalar EffRate12<2, double(*)(const double*, void*)>::
        find_max(std::vector<double> parameters){
    double lnEA = parameters[0], T = parameters[1];
    double EA = std::exp(lnEA);
    double mA = _IS_masses[0];
    double m1 = _FS_masses[0], m2 = _FS_masses[1];
    EA = std::max(std::max(EA, mA), m1+m2);
    double pA = std::sqrt(EA*EA - mA*mA);
    double lnxmin = std::log((.01+m1)/(EA+pA));
    double lnxmax = std::log(1.-m2/(EA+pA));
    double lnsinmin = std::min(std::log(.01/EA), -9.);
    double lnsinmax = 0.;
    auto dR_dlnxdlny = [EA, T, mA, m1, m2, lnxmin, lnxmax, lnsinmin, lnsinmax, this](const double *x){
        if (x[0]<lnxmin||x[0]>lnxmax||x[1]<lnsinmin||x[1]>lnsinmax) return 0.;
        double X[2] = {std::exp(x[0]), std::exp(x[1])};
        double Jacobian = X[0]*X[1];
        double params[5] = {EA, T, mA, m1, m2};
        return this->_f(X, params)*Jacobian;
    };
    auto loc = MC_maximize(dR_dlnxdlny, 2, 
                               {{lnxmin, lnxmax}, {lnsinmin, lnsinmax}}, 
                               2000);
    double xloc[2] = {loc[0],loc[1]};
    return scalar{std::log(1.5*dR_dlnxdlny(xloc))};
}


// Calculate Rate Here! qhat not included yet
template <>
scalar EffRate12<2, double(*)(const double*, void*)>::
        calculate_scalar(std::vector<double> parameters){
    double lnEA = parameters[0], T = parameters[1];
    double EA = std::exp(lnEA);
    double mA = _IS_masses[0];
    double m1 = _FS_masses[0], m2 = _FS_masses[1];
    EA = std::max(std::max(EA, mA), m1+m2);
    double pA = std::sqrt(EA*EA - mA*mA);
    double lnxmin = std::log((.01+m1)/(EA+pA));
    double lnxmax = std::log(1.-m2/(EA+pA));
    double lnsinmin = std::min(std::log(.01/EA), -9.);
    double lnsinmax = 0.;
    auto dR_dlnxdlny = [EA, T, mA, m1, m2, lnxmin, lnxmax, lnsinmin, lnsinmax, this](const double *x){
        if (x[0]<lnxmin||x[0]>lnxmax||x[1]<lnsinmin||x[1]>lnsinmax) return 0.;
        double X[2] = {std::exp(x[0]), std::exp(x[1])};
        double Jacobian = X[0]*X[1];
        double params[5] = {EA, T, mA, m1, m2};
        return this->_f(X, params)*Jacobian;
    };
    double val;

    double xmin[2] = {lnxmin, lnsinmin};
    double xmax[2] = {lnxmax, lnsinmax};
    double err;
    val = vegas(dR_dlnxdlny, 2, xmin, xmax, err, 3000); 
    return scalar{val*_degen};
}


template <size_t N, typename F>
fourvec EffRate12<N, F>::calculate_fourvec(std::vector<double> parameters){
    return fourvec::unity();
}
template <size_t N, typename F>
tensor EffRate12<N, F>::calculate_tensor(std::vector<double> parameters){
    return tensor::unity();
}


template class Rate<HS2PP,2,2,double(*)(const double, void*)>; 
template class Rate<HS2PPP,2,2,double(*)(const double*, void*)>; 
template class EffRate12<2,double(*)(const double*, void*)>; 
