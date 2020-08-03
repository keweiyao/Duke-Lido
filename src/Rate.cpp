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
Rate<HS2HS, 2, 2, double(*)(const double, void *)>::
    Rate(std::string Name, std::string configfile, double(*f)(const double, void *)):
StochasticBase<2>(Name+"/rate", configfile),
X(std::make_shared<Xsection<HS2HS, 2, double(*)(const double, void *)>>(Name, configfile, f) )
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
    _mass = tree.get<double>("mass");
    _degen = tree.get<double>("degeneracy");
    _active = (tree.get<std::string>("<xmlattr>.status")=="active")?true:false;

}

template <>
Rate<HS2QQbar, 2, 2, double(*)(const double, void *)>::
    Rate(std::string Name, std::string configfile, double(*f)(const double, void *)):
StochasticBase<2>(Name+"/rate", configfile),
X(std::make_shared<Xsection<HS2QQbar, 2, double(*)(const double, void *)>>(Name, configfile, f) )
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
    _mass = tree.get<double>("mass");
    _degen = tree.get<double>("degeneracy");
    _active = (tree.get<std::string>("<xmlattr>.status")=="active")?true:false;
}

template <>
Rate<HS2HHS, 2, 2, double(*)(const double*, void *)>::
    Rate(std::string Name, std::string configfile, double(*f)(const double*, void *)):
StochasticBase<2>(Name+"/rate", configfile),
X(std::make_shared<Xsection<HS2HHS, 2, double(*)(const double*, void *)>>(Name, configfile, f) )
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
    _mass = tree.get<double>("mass");
    _degen = tree.get<double>("degeneracy");
    _active = (tree.get<std::string>("<xmlattr>.status")=="active")?true:false;
}

template <>
Rate<HHS2HS, 2, 4, double(*)(const double*, void *)>::
    Rate(std::string Name, std::string configfile, double(*f)(const double*, void *)):
StochasticBase<2>(Name+"/rate", configfile),
X(std::make_shared<Xsection<HHS2HS, 4, double(*)(const double*, void *)>>(Name, configfile, f) )
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
    _mass = tree.get<double>("mass");
    _degen = tree.get<double>("degeneracy");
    _active = (tree.get<std::string>("<xmlattr>.status")=="active")?true:false;
}

/*****************************************************************/
/*************************Sample dR ******************************/
/*****************************************************************/
/*------------------Implementation for 2 -> 2--------------------*/
template <>
void Rate<HS2HS, 2, 2, double(*)(const double, void *)>::
        sample(std::vector<double> parameters,
            std::vector< fourvec > & final_states){
    double lnE = parameters[0], T = parameters[1];
    double E = std::exp(lnE);
    double M2 = _mass*_mass;
    double v1 = std::sqrt(1. - M2/E/E);
    double Q2 = cut*t_channel_mD2->get_mD2(T);
    double sqrtsmin = std::sqrt(M2+.5*Q2 + std::sqrt(M2*Q2+.25*Q2*Q2));
    double lnsqrtsmin = std::log(sqrtsmin);
    auto dR_dxdy = [E, T, lnsqrtsmin, v1, this](const double * x){
        double M = this->_mass;
        double lnsqrts = x[0], costheta = x[1];
        if (lnsqrts<lnsqrtsmin||costheta>=1.||costheta<=-1.) return 0.;
        double s = std::exp(2.*lnsqrts);
        double E2 = (s-M*M)/((1. - v1*costheta)*2.*E);
        double Jacobian = 2.*E2;
        double Xtot = this->X->GetZeroM({lnsqrts,T}).s;
        return Jacobian/E*E2*std::exp(-E2/T)*(s-M*M)*2*Xtot/16./M_PI/M_PI;
    };
    bool status = true;
    double fmax = std::exp(StochasticBase<2>::GetFmax(parameters).s);
    auto res = sample_nd(dR_dxdy, 2, 
                      {{lnsqrtsmin, std::log(sqrtsmin+100*E*T)}, {-1., 1.}},
                      fmax, status);
    double lnsqrts = res[0], costheta = res[1];
    double s = std::exp(2.*lnsqrts);
    double E2 = (s-_mass*_mass)/((1. - v1*costheta)*2.*E);
    double sintheta = std::sqrt(1. - costheta*costheta);

    double tmin = -std::pow(s-_mass*_mass, 2)/s, 
           tmax = -cut*t_channel_mD2->get_mD2(T);
    X->sample({lnsqrts, T}, final_states);
    // give incoming partilce a random phi angle
    double phi = Srandom::dist_phi(Srandom::gen);
    // com velocity
    double vcom[3] = { E2*sintheta/(E2+E)*std::cos(phi),
               E2*sintheta/(E2+E)*std::sin(phi),
               (E2*costheta+v1*E)/(E2+E)    };
    /*
    FS now is in Z-oriented CoM frame
    1) FS.rotate_back
    2) FS.boost_back
    */
    fourvec p1{E, 0, 0, v1*E};
    auto p1com = p1.boost_to(vcom[0], vcom[1], vcom[2]);
    for(auto & p: final_states){
        p = p.rotate_back(p1com);
        p = p.boost_back(vcom[0], vcom[1], vcom[2]);
    }
}


template <>
void Rate<HS2QQbar, 2, 2, double(*)(const double, void *)>::
        sample(std::vector<double> parameters,
            std::vector< fourvec > & final_states){
    double lnE = parameters[0], T = parameters[1];
    double E = std::exp(lnE);
    auto dR_dlnsqrts_dcostheta = [E, T, this](const double * x){
        double lnsqrts = x[0], costheta = x[1];
        double sqrts = std::exp(lnsqrts);
        double s = sqrts*sqrts;
        double E2 = s/(2.*E*(1. - costheta));
        double Jacobian = 2*E2;
        double Xtot = this->X->GetZeroM({lnsqrts,T}).s;
        return Jacobian/E*E2*std::exp(-E2/T)*s*2*Xtot/16./M_PI/M_PI;
    };
    double xmin[2] = {std::log(2.*_mass), -1.};
    double xmax[2] = {std::log(2.*_mass*100), 1.};

    bool status = true;
    auto res = sample_nd(dR_dlnsqrts_dcostheta, 2, 
               {{std::log(2.*_mass), std::log(100.*2.*_mass)}, 
                {-1., 1.}},
               std::exp(StochasticBase<2>::GetFmax(parameters).s), 
               status);

    double lnsqrts = res[0], costheta = res[1];
    double E2 = std::exp(2.*lnsqrts)/(2.*E*(1. - costheta));
    double sintheta = std::sqrt(1. - costheta*costheta);
    X->sample({lnsqrts, T}, final_states);

    // give incoming partilce a random phi angle
    double phi = Srandom::dist_phi(Srandom::gen);
    // com velocity
    double vcom[3] = { E2*sintheta/(E2+E)*std::cos(phi),
               E2*sintheta/(E2+E)*std::sin(phi),
              (E2*costheta+E)/(E2+E)    };
    /*
    FS now is in Z-oriented CoM frame
    1) FS.rotate_back
    2) FS.boost_back
    */
    fourvec p1{E, 0, 0, E};
    auto p1com = p1.boost_to(vcom[0], vcom[1], vcom[2]);
    for(auto & p: final_states){
        p = p.rotate_back(p1com);
        p = p.boost_back(vcom[0], vcom[1], vcom[2]);
    }
}

/*------------------Implementation for 2 -> 3--------------------*/
template <>
void Rate<HS2HHS, 2, 2, double(*)(const double*, void *)>::
        sample(std::vector<double> parameters,
            std::vector< fourvec > & final_states){
    double lnE = parameters[0], T = parameters[1];
    double E = std::exp(lnE);
    double M2 = _mass*_mass;
    double v1 = std::sqrt(1. - M2/E/E);
    double Q2 = cut*t_channel_mD2->get_mD2(T);
    double sqrtsmin = std::sqrt(M2+.5*Q2 + std::sqrt(M2*Q2+.25*Q2*Q2));
    double lnsqrtsmin = std::log(sqrtsmin);
    auto dR_dxdy = [E, T, lnsqrtsmin, v1, this](const double * x){
        double M = this->_mass;
        double lnsqrts = x[0], costheta = x[1];
        if (lnsqrts<lnsqrtsmin||costheta>=1.||costheta<=-1.) return 0.;
        double s = std::exp(2.*lnsqrts);
        double E2 = (s-M*M)/((1. - v1*costheta)*2.*E);
        double Jacobian = 2.*E2;
        double Xtot = this->X->GetZeroM({lnsqrts,T}).s;
        return Jacobian/E*E2*std::exp(-E2/T)*(s-M*M)*2*Xtot/16./M_PI/M_PI;
    };
    bool status = true;
    auto res = sample_nd(dR_dxdy, 2, 
               {{std::log(sqrtsmin), std::log(sqrtsmin*10)}, 
               {-1., 1.}}, 
               std::exp(StochasticBase<2>::GetFmax(parameters).s), status);
    double lnsqrts = res[0], costheta = res[1];
    double sintheta = std::sqrt(1. - costheta*costheta);
    double phi = Srandom::dist_phi(Srandom::gen);
    double cosphi = std::cos(phi), sinphi = std::sin(phi);
    double s = std::exp(2.*lnsqrts);
    double E2 = (s-M2)/((1. - v1*costheta)*2.*E);
    double vcom[3] = { E2*sintheta/(E2+E)*cosphi,
                        E2*sintheta/(E2+E)*sinphi,
                        (E2*costheta+v1*E)/(E2+E) };
    X->sample({lnsqrts, T}, final_states);

    fourvec p1{E, 0, 0, v1*E};
    auto p1com = p1.boost_to(vcom[0], vcom[1], vcom[2]);
    for(auto & p: final_states){
        p = p.rotate_back(p1com);
        p = p.boost_back(vcom[0], vcom[1], vcom[2]);
    }
}

/*------------------Implementation for 3 -> 2--------------------*/
template <>
void Rate<HHS2HS, 2, 4, double(*)(const double*, void *)>::
        sample(std::vector<double> parameters,
            std::vector< fourvec > & final_states){
    double lnE = parameters[0], T = parameters[1];
    double E = std::exp(lnE);
    double M = _mass;
    double M2 = M*M;
    // sample dR
    // x are: k, E2, cosk, cos2, phi2
    auto code = [E, T, this](const double * x){
        double M = this->_mass;
        double M2 = M*M;
        double k = x[0], E2 = x[1], cosk = x[2], cos2 = x[3], phi2 = x[4];
        double sink = std::sqrt(1.-cosk*cosk), sin2 = std::sqrt(1.-cos2*cos2);
        double cosphi2 = std::cos(phi2), sinphi2 = std::sin(phi2);
        double v1 = std::sqrt(1. - M*M/E/E);
        fourvec p1mu{E, 0, 0, v1*E};
        fourvec p2mu{E2, E2*sin2*cosphi2, E2*sin2*sinphi2, E2*cos2};
        fourvec kmu{k, k*sink, 0., k*cosk};
        fourvec Ptot = p1mu+p2mu+kmu, P12 = p1mu+p2mu, P1k = p1mu+kmu;
        double s = dot(Ptot, Ptot), s12 = dot(P12, P12), s1k = dot(P1k, P1k);
        double sqrts = std::sqrt(s), sqrts12 = std::sqrt(s12), sqrts1k = std::sqrt(s1k);
        double lnsqrts = std::log(sqrts);
        double v12[3] = { P12.x()/P12.t(), P12.y()/P12.t(), P12.z()/P12.t() };
        double xinel = (s12-M2)/(s-M2), yinel = (s1k/s-M2/s12)/(1.-s12/s)/(1.-M2/s12);
        // interp Xsection
        double Xtot = std::exp(this->X->GetZeroM({lnsqrts, T, xinel, yinel}).s);
        return std::exp(-(k+E2)/T)*k*E2*Xtot/E/8./std::pow(2.*M_PI, 5);
    };
    bool status = true;
    auto res = sample_nd(code, 5, {{0.0*T, 10.0*T}, {0.0*T, 10.0*T}, {-1., 1.}, {-1., 1.}, {0., 2.*M_PI}}, std::exp(StochasticBase<2>::GetFmax(parameters).s), status);
    // sample Xsection
    double k = res[0];
    double E2 = res[1];
    double cosk = res[2];
    double cos2 = res[3];
    double sink = std::sqrt(1.-cosk*cosk), sin2 = std::sqrt(1.-cos2*cos2);
    double phik = Srandom::dist_phi(Srandom::gen);// randomed phi_k
    double phi2 = res[4]+phik; // phi2 is relative to phi_k
    double cosphi2 = std::cos(phi2), sinphi2 = std::sin(phi2);
    double v1 = std::sqrt(1. - M*M/E/E);
    fourvec p1mu{E, 0, 0, v1*E};
    fourvec p2mu{E2, E2*sin2*cosphi2, E2*sin2*sinphi2, E2*cos2};
    fourvec kmu{k, k*sink*std::cos(phik), k*sink*std::sin(phik), k*cosk};
    fourvec Ptot = p1mu+p2mu+kmu, P12 = p1mu+p2mu, P1k = p1mu+kmu;
    double s = dot(Ptot, Ptot), s12 = dot(P12, P12), s1k = dot(P1k, P1k);
    double sqrts = std::sqrt(s), sqrts12 = std::sqrt(s12), sqrts1k = std::sqrt(s1k);
    double lnsqrts = std::log(sqrts);
    double v12[3] = { P12.x()/P12.t(), P12.y()/P12.t(), P12.z()/P12.t() };
    double xinel = (s12-M2)/(s-M2), yinel = (s1k/s-M2/s12)/(1.-s12/s)/(1.-M2/s12);
    X->sample({lnsqrts,T,xinel,yinel}, final_states);

    // p1 and k in 12-com frame
    auto p1_in12 = p1mu.boost_to(v12[0], v12[1], v12[2]);
    auto k_in12 = kmu.boost_to(v12[0], v12[1], v12[2]);
    // |p1|
    double p1abs = std::sqrt(p1_in12.x()*p1_in12.x() + p1_in12.y()*p1_in12.y() + p1_in12.z()*p1_in12.z());
    // X-frame z-direction
    double zdir[3] = {p1_in12.x()/p1abs, p1_in12.y()/p1abs, p1_in12.z()/p1abs};
    // project k onto z-dir
    double kdotz = zdir[0]*k_in12.x() + zdir[1]*k_in12.y() + zdir[2]*k_in12.z();
    // get kperp_to-zdir
    double kperp[3] = {k_in12.x() - kdotz*zdir[0], k_in12.y() - kdotz*zdir[1], k_in12.z() - kdotz*zdir[2]};
    double kperpabs = std::sqrt(kperp[0]*kperp[0]+kperp[1]*kperp[1]+kperp[2]*kperp[2]);
    // define x-dir
    double xdir[3] = {kperp[0]/kperpabs, kperp[1]/kperpabs, kperp[2]/kperpabs};
    // y-dir = z-dir \cross x-dir
    double ydir[3] = {zdir[1]*xdir[2]-zdir[2]*xdir[1], xdir[0]*zdir[2]-xdir[2]*zdir[0], zdir[0]*xdir[1]-xdir[0]*zdir[1]};
    for (auto & p : final_states){
        // rotate it back
        double px = p.x()*xdir[0] + p.y()*ydir[0] + p.z()*zdir[0];
        double py = p.x()*xdir[1] + p.y()*ydir[1] + p.z()*zdir[1];
        double pz = p.x()*xdir[2] + p.y()*ydir[2] + p.z()*zdir[2];
        p.a[1] = px;
        p.a[2] = py;
        p.a[3] = pz;
        p = p.boost_back(v12[0], v12[1], v12[2]);
    }
    final_states.push_back(kmu);
}

/*****************************************************************/
/*************************Find dR_max ****************************/
/*****************************************************************/
/*------------------Implementation for 2 -> 2--------------------*/
template <>
scalar Rate<HS2HS, 2, 2, double(*)(const double, void*)>::
        find_max(std::vector<double> parameters){
    double lnE = parameters[0], T = parameters[1];
    double E = std::exp(lnE);
    double Q2 = cut*t_channel_mD2->get_mD2(T);
    double M2 = _mass*_mass;
    double sqrtsmin = std::sqrt(M2+.5*Q2 + std::sqrt(M2*Q2+.25*Q2*Q2));
    double lnsqrtsmin = std::log(sqrtsmin);
    auto minus_dR_dxdy = [E, T, lnsqrtsmin, this](const double * x){
        double M = this->_mass;
        double v1 = std::sqrt(1. - M*M/E/E);
        double lnsqrts = x[0], costheta = x[1];
        if (lnsqrts<lnsqrtsmin||costheta>=1.||costheta<=-1.) return 0.;
        double s = std::exp(2.*lnsqrts);
        double E2 = (s-M*M)/((1. - v1*costheta)*2.*E);
        double Jacobian = 2.*E2;
        double Xtot = this->X->GetZeroM({lnsqrts,T}).s;
        return -Jacobian/E*E2*std::exp(-E2/T)*(s-M*M)*2*Xtot/16./M_PI/M_PI;
    };
    auto val = -minimize_nd(minus_dR_dxdy, 2, 
                            {lnsqrtsmin*2, 0.}, {lnsqrtsmin*.5, 0.2}, 1000, 1e-8);
    return scalar{std::log(1.5*val)};
}

template <>
scalar Rate<HS2QQbar, 2, 2, double(*)(const double, void*)>::
        find_max(std::vector<double> parameters){
    double lnE = parameters[0], T = parameters[1];
    double E = std::exp(lnE);
    auto dR_dxdy = [E, T, this](const double * x){
        double lnsqrts = x[0], costheta = x[1];
        double sqrts = std::exp(lnsqrts);
        if (sqrts<2.*this->_mass || costheta >= 1. || costheta <= -1.) 
            return 0.;
        double s = sqrts*sqrts;
        double E2 = s/(2.*E*(1. - costheta));
        double Jacobian = 2*E2;
        double Xtot = this->X->GetZeroM({lnsqrts,T}).s;
        return -Jacobian/E*E2*std::exp(-E2/T)*s*2*Xtot/16./M_PI/M_PI;
    };
    auto val = -minimize_nd(dR_dxdy, 2, {std::log(3.*_mass), 0.}, 
                            {0.2, 0.2}, 1000, 1e-8)*1.5;
    return scalar{std::log(val)};
}
/*------------------Implementation for 2 -> 3--------------------*/
template <>
scalar Rate<HS2HHS, 2, 2, double(*)(const double*, void*)>::
        find_max(std::vector<double> parameters){
    double lnE = parameters[0], T = parameters[1];
    double E = std::exp(lnE);
    double M2 = _mass*_mass;
    double v1 = std::sqrt(1. - M2/E/E);
    double Q2 = cut*t_channel_mD2->get_mD2(T);
    double sqrtsmin = std::sqrt(M2+.5*Q2 + std::sqrt(M2*Q2+.25*Q2*Q2));
    double lnsqrtsmin = std::log(sqrtsmin);
    auto dR_dxdy = [E, T, lnsqrtsmin, v1, this](const double * x){
        double M = this->_mass;
        double lnsqrts = x[0], costheta = x[1];
        if (lnsqrts<lnsqrtsmin||costheta>=1.||costheta<=-1.) return 0.;
        double s = std::exp(2.*lnsqrts);
        double E2 = (s-M*M)/((1. - v1*costheta)*2.*E);
        double Jacobian = 2.*E2;
        double Xtot = this->X->GetZeroM({lnsqrts,T}).s;
        return -Jacobian/E*E2*std::exp(-E2/T)*(s-M*M)*2*Xtot/16./M_PI/M_PI;
    };
    auto val = -minimize_nd(dR_dxdy, 2, {lnsqrtsmin*2, 0.}, {lnsqrtsmin*.5, 0.2}, 1000, 1e-8)*3;
    return scalar{std::log(val)};
}
/*------------------Implementation for 3 -> 2--------------------*/
template <>
scalar Rate<HHS2HS, 2, 4, double(*)(const double*, void*)>::
        find_max(std::vector<double> parameters){
    double lnE = parameters[0], T = parameters[1];
    double E = std::exp(lnE);
    // x are: k, E2, cosk, cos2, phi2
    auto code = [E, T, this](const double * x){
        double M = this->_mass;
        double M2 = M*M;
        double k = x[0], E2 = x[1], cosk = x[2], cos2 = x[3], phi2 = x[4];
        if (k<0.01*T || E2<0.01*T || k>8*T || E2>8*T || std::abs(cosk)>.999|| std::abs(cos2)>.999|| phi2<0. || phi2 > 2.*M_PI) return 0.;
        double sink = std::sqrt(1.-cosk*cosk), sin2 = std::sqrt(1.-cos2*cos2);
        double cosphi2 = std::cos(phi2), sinphi2 = std::sin(phi2);
        double v1 = std::sqrt(1. - M*M/E/E);
        fourvec p1mu{E, 0, 0, v1*E};
        fourvec p2mu{E2, E2*sin2*cosphi2, E2*sin2*sinphi2, E2*cos2};
        fourvec kmu{k, k*sink, 0., k*cosk};
        fourvec Ptot = p1mu+p2mu+kmu, P12 = p1mu+p2mu, P1k = p1mu+kmu;
        double s = dot(Ptot, Ptot), s12 = dot(P12, P12), s1k = dot(P1k, P1k);
        double sqrts = std::sqrt(s), sqrts12 = std::sqrt(s12), sqrts1k = std::sqrt(s1k);
        double lnsqrts = std::log(sqrts);
        double v12[3] = { P12.x()/P12.t(), P12.y()/P12.t(), P12.z()/P12.t() };
        double xinel = (s12-M2)/(s-M2), yinel = (s1k/s-M2/s12)/(1.-s12/s)/(1.-M2/s12);
        // interp Xsection
        double Xtot = std::exp(this->X->GetZeroM({lnsqrts, T, xinel, yinel}).s);
        return std::exp(-(k+E2)/T)*k*E2*Xtot/E/8./std::pow(2.*M_PI, 5);
    };
    auto loc = MC_maximize(code, 5,
                {{T, T*2}, {T, T*2},{-.5, .5}, {-0.5, 0.5}, {1.5,2.}}, 500);
    double xloc[5] = {loc[0],loc[1],loc[2],loc[3],loc[4]};
    auto val = code(xloc)*2.;
    return scalar{std::log(val)};
}


/*****************************************************************/
/*************************Integrate dR ***************************/
/*****************************************************************/
/*------------------Implementation for 2 -> 2--------------------*/
template <>
scalar Rate<HS2HS, 2, 2, double(*)(const double, void*)>::
        calculate_scalar(std::vector<double> parameters){
    double lnE = parameters[0], T = parameters[1];
    double E = std::exp(lnE);
    double M2 = _mass*_mass;
    double v1 = std::sqrt(1. - M2/E/E);
    double Q2 = cut*t_channel_mD2->get_mD2(T);
    double sqrtsmin = std::sqrt(M2+.5*Q2 + std::sqrt(M2*Q2+.25*Q2*Q2));
    auto code = [E, T, v1, this](const double * x){
        double M = this->_mass;
        double lnsqrts = x[0], costheta = x[1];
        double s = std::exp(2.*lnsqrts);
        double E2 = (s-M*M)/((1. - v1*costheta)*2.*E);
        double Jacobian = 2.*E2;
        std::vector<double> res{0.};
        double Xtot = this->X->GetZeroM({lnsqrts,T}).s;
        res[0] = Jacobian/E*E2*std::exp(-E2/T)*(s-M*M)*2*Xtot/16./M_PI/M_PI;
        return res;
    };
    double xmin[2] = {std::log(sqrtsmin), -1.};
    double xmax[2] = {std::log(sqrtsmin+100*E*T),1.};
    double err;
    auto val = quad_nd(code, 2, 1, xmin, xmax, err);
    return scalar{_degen*val[0]};
}

template <>
scalar Rate<HS2QQbar, 2, 2, double(*)(const double, void*)>::
        calculate_scalar(std::vector<double> parameters){
    double lnE = parameters[0], T = parameters[1];
    double E = std::exp(lnE);
    auto code = [E, T, this](const double * x){
        double lnsqrts = x[0], costheta = x[1];
        double sqrts = std::exp(lnsqrts);
        double s = sqrts*sqrts;
        double E2 = s/(2.*E*(1. - costheta));
        double Jacobian = 2*E2;
        double Xtot = this->X->GetZeroM({lnsqrts,T}).s;
        std::vector<double> res{Jacobian/E*E2*std::exp(-E2/T)*s*2*Xtot/16./M_PI/M_PI};
        return res;
    };
    double xmin[2] = {std::log(2.*_mass), -1.};
    double xmax[2] = {std::log(2.*_mass*100), 1.};
    double err;
    auto val = quad_nd(code, 2, 1, xmin, xmax, err);
    return scalar{_degen*val[0]};
}
/*------------------Implementation for 2 -> 3--------------------*/
template <>
scalar Rate<HS2HHS, 2, 2, double(*)(const double*, void*)>::
        calculate_scalar(std::vector<double> parameters){
    double lnE = parameters[0], T = parameters[1];
    double E = std::exp(lnE);
    double M2 = _mass*_mass;
    double v1 = std::sqrt(1. - M2/E/E);
    double Q2 = cut*t_channel_mD2->get_mD2(T);
    double sqrtsmin = std::sqrt(M2+.5*Q2 + std::sqrt(M2*Q2+.25*Q2*Q2));
    double lnsqrtsmin = std::log(sqrtsmin);
    auto dR_dxdy = [E, T, lnsqrtsmin, v1, this](const double * x){
        double M = this->_mass;
        double lnsqrts = x[0], costheta = x[1];
        double s = std::exp(2.*lnsqrts);
        double E2 = (s-M*M)/((1. - v1*costheta)*2.*E);
        double Jacobian = 2.*E2;
        double Xtot = this->X->GetZeroM({lnsqrts,T}).s;
        std::vector<double> res{Jacobian/E*E2*std::exp(-E2/T)*(s-M*M)*2*Xtot/16./M_PI/M_PI};
        return res;
    };
    double xmin[2] = {std::log(sqrtsmin), -1.};
    double xmax[2] = {std::log(sqrtsmin*10),1.};
    double err;
    auto val = quad_nd(dR_dxdy, 2, 1, xmin, xmax, err);
    return scalar{_degen*val[0]};
}

/*------------------Implementation for 3 -> 2--------------------*/
template <>
scalar Rate<HHS2HS, 2, 4, double(*)(const double*, void*)>::
        calculate_scalar(std::vector<double> parameters){
    double lnE = parameters[0], T = parameters[1];
    double E = std::exp(lnE);
    // x are: k, E2, cosk, cos2, phi2
    auto code = [E, T, this](const double * x){
        double M = this->_mass;
        double M2 = M*M;
        double k = x[0], E2 = x[1], cosk = x[2], cos2 = x[3], phi2 = x[4];
        double sink = std::sqrt(1.-cosk*cosk), sin2 = std::sqrt(1.-cos2*cos2);
        double cosphi2 = std::cos(phi2), sinphi2 = std::sin(phi2);
        double v1 = std::sqrt(1. - M*M/E/E);
        fourvec p1mu{E, 0, 0, v1*E};
        fourvec p2mu{E2, E2*sin2*cosphi2, E2*sin2*sinphi2, E2*cos2};
        fourvec kmu{k, k*sink, 0., k*cosk};
        fourvec Ptot = p1mu+p2mu+kmu, P12 = p1mu+p2mu, P1k = p1mu+kmu;
        double s = dot(Ptot, Ptot), s12 = dot(P12, P12), s1k = dot(P1k, P1k);
        double sqrts = std::sqrt(s), sqrts12 = std::sqrt(s12), sqrts1k = std::sqrt(s1k);
        double lnsqrts = std::log(sqrts);
        double v12[3] = { P12.x()/P12.t(), P12.y()/P12.t(), P12.z()/P12.t() };
        double xinel = (s12-M2)/(s-M2), yinel = (s1k/s-M2/s12)/(1.-s12/s)/(1.-M2/s12);
        // interp Xsection
        double Xtot = std::exp(this->X->GetZeroM({lnsqrts, T, xinel, yinel}).s);
        return std::exp(-(k+E2)/T)*k*E2*Xtot/E/8./std::pow(2.*M_PI, 5);
    };
    double xmin[5] = {0.,   0.,  -1., -1, 0.};
    double xmax[5] = {10*T, 10*T, 1., 1., 2.*M_PI};
    double error;
    double res = vegas(code, 5, xmin, xmax, error);
    return scalar{_degen*res};
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
fourvec Rate<HS2HS, 2, 2, double(*)(const double, void*)>::
        calculate_fourvec(std::vector<double> parameters){
    double lnE = parameters[0], T = parameters[1];
    double E = std::exp(lnE);
    double Q2 = cut*t_channel_mD2->get_mD2(T);
    double M2 = _mass*_mass;
    double v1 = std::sqrt(1. - M2/E/E);
    double sqrtsmin = std::sqrt(M2+.5*Q2 + std::sqrt(M2*Q2+.25*Q2*Q2));
    auto code = [E, T, v1,this](const double * x){
        double M = this->_mass;
        double lnsqrts = x[0], costheta = x[1];
        double sintheta = std::sqrt(1.-costheta*costheta);
        double s = std::exp(2.*lnsqrts);
        double E2 = (s-M*M)/((1. - v1*costheta)*2.*E);
        double Jacobian = 2.*E2;
        double vcom[3] = {E2*sintheta/(E2+E), 0., (E2*costheta+v1*E)/(E2+E)};
        fourvec p1{E, 0, 0, v1*E};
        // A vector in p1z(com)-oriented com frame
        auto fmu0 = this->X->GetFirstM({lnsqrts,T});
        // rotate it back from p1z(com)-oriented com frame
        auto fmu1 = fmu0.rotate_back(p1.boost_to(vcom[0], vcom[1], vcom[2]));
        // boost back to the matter frame
        auto fmu2 = fmu1.boost_back(vcom[0], vcom[1], vcom[2]);
        double common = Jacobian/E*E2*std::exp(-E2/T)*(s-M*M)*2./16./M_PI/M_PI;
        fmu2 = fmu2 * common;
        // Set tranverse to zero due to azimuthal symmetry;
        std::vector<double> res{fmu2.t(), fmu2.z()};
        return res;
    };
    double xmin[2] = {0., -1.};
    double xmax[2] = {10.*T, 1.};
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
tensor Rate<HS2HS, 2, 2, double(*)(const double, void*)>::
        calculate_tensor(std::vector<double> parameters){
    double lnE = parameters[0], T = parameters[1];
    double E = std::exp(lnE);
    double M2 = _mass*_mass;
    double v1 = std::sqrt(1. - M2/E/E);
    double Q2 = cut*t_channel_mD2->get_mD2(T);
    double sqrtsmin = std::sqrt(M2+.5*Q2 + std::sqrt(M2*Q2+.25*Q2*Q2));
    auto code = [E, T, v1,this](const double * x){
        double M = this->_mass;
        double lnsqrts = x[0], costheta = x[1], phi=x[2];
        double s = std::exp(2.*lnsqrts);
        double E2 = (s-M*M)/((1. - v1*costheta)*2.*E);
        double Jacobian = 2.*E2;
        double sintheta = std::sqrt(1. - costheta*costheta);
        double vcom[3] = {E2*sintheta/(E2+E)*std::cos(phi), 
                    E2*sintheta/(E2+E)*std::sin(phi),
                 (E2*costheta+v1*E)/(E2+E)};
        fourvec p1{E, 0, 0, v1*E};
        // A vector in p1z(com)-oriented com frame
        auto fmunu0 = this->X->GetSecondM({lnsqrts,T});
        // rotate it back from p1z(com)-oriented com frame
        auto fmunu1 = fmunu0.rotate_back(p1.boost_to(vcom[0], vcom[1], vcom[2]));
        // boost back to the matter frame
        auto fmunu2 = fmunu1.boost_back(vcom[0], vcom[1], vcom[2]);
        double common = Jacobian/E*E2*std::exp(-E2/T)*(s-M*M)*2./32./std::pow(M_PI, 3);
        fmunu2 = fmunu2 * common;
        // Set tranverse to zero due to azimuthal symmetry;
        std::vector<double> res{fmunu2.T[0][0], fmunu2.T[1][1],
                                fmunu2.T[2][2], fmunu2.T[3][3]};
        return res;
    };
    double xmin[3] = {std::log(sqrtsmin), -1., -M_PI};
    double xmax[3] = {std::log(sqrtsmin+6*E*T), 1., M_PI};
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
   _mass = tree.get<double>("mass");
   _active = (tree.get<std::string>("<xmlattr>.status")=="active")?true:false;
}

// Sample Final states
template<>
void EffRate12<2, double(*)(const double*, void *)>::
    sample(std::vector<double> parameters, std::vector< fourvec > & final_states){
    double lnE = parameters[0], T = parameters[1];
    double E = std::exp(lnE);
    double pabs = std::sqrt(E*E - _mass*_mass);

    auto dR_dlnxdlny = [E, T, this](const double *x){
        double X[2] = {std::exp(x[0]), x[1]};
        double Jacobian = X[0];
        double params[3] = {E, T, this->_mass};
        return this->_f(X, params)*Jacobian;
    };
    bool status=true;
    double fmax = StochasticBase<2>::GetFmax(parameters).s;
    auto res = sample_nd(dR_dlnxdlny, 2, 
                        {{-9. ,0.}, {0., 1.}}, 
                         fmax,
                         status);

    double x = std::exp(res[0]), y = res[1];
    double k0 = x*pabs;
    double kT = x*y*E;
    double phi = Srandom::dist_phi(Srandom::gen);
    double kx = kT*std::cos(phi), ky = kT*std::sin(phi);
    double kz = std::sqrt(k0*k0-kT*kT);
    double pznew = pabs-kz;
    double Enew = std::sqrt(pznew*pznew + kT*kT + _mass*_mass);

    final_states.resize(2);
    final_states[0] = fourvec{Enew, -kx, -ky, pznew};
    final_states[1] = fourvec{k0, kx, ky, kz};
}

// Find the max of dR/dx/dy
template <>
scalar EffRate12<2, double(*)(const double*, void*)>::
        find_max(std::vector<double> parameters){
    double lnE = parameters[0], T = parameters[1];
    double E = std::exp(lnE);
    double pabs = std::sqrt(E*E-_mass*_mass);
    auto dR_dlnxdlny = [E, T, this](const double * x){
        double X[2] = {std::exp(x[0]), x[1]};
        if (X[0]<=0.||X[0]>=1.||X[1]<=0.||X[1]>=1.) return 0.;
        double Jacobian = X[0];
        double params[3] = {E, T, this->_mass};
        return this->_f(X, params)*Jacobian;
    };
    auto loc = MC_maximize(dR_dlnxdlny, 2, 
                               {{-9.,0.}, {0,1}}, 
                               2000);
    double xloc[2] = {loc[0],loc[1]};
    return scalar{dR_dlnxdlny(xloc)*1.1};
}


// Calculate Rate Here! qhat not included yet
template <>
scalar EffRate12<2, double(*)(const double*, void*)>::
        calculate_scalar(std::vector<double> parameters){
    double lnE = parameters[0], T = parameters[1];
    double E = std::exp(lnE);
    double pabs = std::sqrt(E*E-_mass*_mass);
    // dR/dx/dy to be intergated
    auto dR_dlnxdlny = [E, T, this](const double * x){
        double X[2] = {std::exp(x[0]), x[1]};
        if (X[0]<=0.||X[0]>=1.||X[1]<=0.||X[1]>=1.) return 0.;
        double Jacobian = X[0];
        double params[3] = {E, T, this->_mass};
        return this->_f(X, params)*Jacobian;
    };
    double val;

    double xmin[2] = {-9., 0.};
    double xmax[2] = {0., 1.};
    double err;
    val = vegas(dR_dlnxdlny, 2, xmin, xmax, err, 3000); 
    return scalar{val};
}


template <size_t N, typename F>
fourvec EffRate12<N, F>::calculate_fourvec(std::vector<double> parameters){
    return fourvec::unity();
}
template <size_t N, typename F>
tensor EffRate12<N, F>::calculate_tensor(std::vector<double> parameters){
    return tensor::unity();
}


//EffRate21 constructor
template<>
EffRate21<2, double(*)(const double*, void *)>::
    EffRate21(std::string Name, std::string configfile, double(*f)(const double *, void *)):
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
    _mass = tree.get<double>("mass");
    _active = (tree.get<std::string>("<xmlattr>.status") == "active") ? true: false;
}

// Sample Final status
template<>
void EffRate21<2, double(*)(const double*, void *)>::
    sample(std::vector<double> parameters, std::vector< fourvec > & final_states){
    double lnE = parameters[0], T = parameters[1];
    double E = std::exp(lnE);
    double pabs = std::sqrt(E*E - _mass*_mass);

    auto dR_dxdy = [E, T, this](const double *x){
        double params[3] = {E, T, this->_mass};
        double result = this->_f(x, params);
        return result;
    };

    bool status=true;
    auto res = sample_nd(dR_dxdy, 2, {{-.99999, .99999}, {0., 1.}}, StochasticBase<2>::GetFmax(parameters).s, status);

    double xp = res[0], y = res[1];
    double x = std::atanh(xp);
    double kz = x*pabs;
    double k0 = std::abs(kz)/std::sqrt(1.01-y*y);
    double kT = y*k0;
    double phi = Srandom::dist_phi(Srandom::gen);
    double kx = kT*std::cos(phi), ky = kT*std::sin(phi);
    double pznew = pabs+kz;
    double Enew = std::sqrt(pznew*pznew + kT*kT + _mass*_mass);

    final_states.resize(2);
    final_states[0] = fourvec{Enew, kx, ky, pznew};
    final_states[1] = fourvec{k0, kx, ky, kz};

    return ;
}


// Find the max of dR/dx/dy
template <>
scalar EffRate21<2, double(*)(const double*, void*)>::
        find_max(std::vector<double> parameters){
    double lnE = parameters[0], T = parameters[1];
    double E = std::exp(lnE);
    auto dR_dxdy = [E, T, this](const double * x){
        // if out of physics range, please return 0.
        if (x[0] >= .99999 || x[0] <= 0. || x[1] <=0 || x[1] >= 1)
            return 0.;
        double params[4] = {E, T, this->_mass};
        double result = this->_f(x, params);
        return result;
    };
    auto loc = MC_maximize(dR_dxdy, 2, {{0,1}, {0,1}}, 1000);
    double xloc[2] = {loc[0],loc[1]};
    
    return scalar{dR_dxdy(xloc)*2};
}




// Calculate the Rate here! qhat not included yet
template<>
scalar EffRate21<2, double(*)(const double*, void*)>::
        calculate_scalar(std::vector<double> parameters){
    double lnE = parameters[0], T = parameters[1];
    double E = std::exp(lnE);

    auto dR_dxdy = [E, T, this](const double *x){
        double params[3] = {E, T, this->_mass};
        double result = this->_f(x, params);
        std::vector<double> res{result};
        return result;
    };

    double xmin[2] = {-.99999, 0.};
    double xmax[2] = {.99999, 1.};
    double err;

    auto val = vegas(dR_dxdy, 2, xmin, xmax, err, 10000);
    return scalar{val};
}

template<size_t N, typename F>
fourvec EffRate21<N, F>::calculate_fourvec(std::vector<double> parameters){
    return fourvec::unity();
}
template<size_t N, typename F>
tensor EffRate21<N, F>::calculate_tensor(std::vector<double> parameters){
    return tensor::unity();    
}


// For 2->2 with leading order pQCD
template class Rate<HS2HS,2,2,double(*)(const double, void*)>; 
template class Rate<HS2QQbar,2,2,double(*)(const double, void*)>; 
template class Rate<HS2HHS,2,2,double(*)(const double*, void*)>; 
template class Rate<HHS2HS,2,4,double(*)(const double*, void*)>;
template class EffRate12<2,double(*)(const double*, void*)>; // 1->2 
template class EffRate21<2,double(*)(const double*, void*)>; // 2->1 
