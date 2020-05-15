#ifndef JET_FINDING_H
#define JET_FINDING_H

#include <fstream>
#include <vector>
#include <string>
#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/FastJet3.h"
#include "lorentz.h"
#include "workflow.h"
#include "TableBase.h"

class MediumResponse{
private:
    // A table for response dot fucntion G^mu,
    // the response is G \cdot Delta P
    std::shared_ptr<TableBase<fourvec, 4>> Gmu;
    std::string name;
    void compute(int start, int end);
public:
    MediumResponse(std::string fname);
    double get_dpT_dydphi(double y, double phi, fourvec Pmu, double vperp, double pTmin);
    void init(std::string);
    void load(std::string);
};

struct current{
    fourvec p;
    double chetas, shetas;
    double cs;
};

struct HadronizeCurrent{
    fourvec G;
    double UdotG, phiu, gamma_perp, v_perp, etas;
};

template <typename T>
std::string to_string(T value)
{
    std::ostringstream os;
    os << value;
    return os.str();
}

template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 2)
{
    std::ostringstream out;
    out << std::setprecision(n) << a_value;
    return out.str();
}

class MyInfo: public fastjet::PseudoJet::UserInfoBase {
    public:
    MyInfo(int pid, int origin, double weight) : _pdg_id(pid), _mc_origin(origin), _w(weight) {};
    
    int pdg_id() const {return _pdg_id;}
    int mc_origin() const {return _mc_origin;}
    double w() const {return _w;}
    private:
    int _pdg_id, _mc_origin;
    double _w;
};
struct Fjet{
    fourvec pmu;
    double R, pT, M2, phi, eta;
    std::vector<double> shape;
    int flavor;
    double sigma;
};
void TestSource(
             MediumResponse MR,
             std::vector<current> jlist,
             std::vector<HadronizeCurrent> slist,
             std::string fname
     );

class JetFinder{
public:
    JetFinder(int Neta, int Nphi, double etamax, bool read_table, std::string stable_path);
    JetFinder(int Neta, int Nphi, double etamax);
    void MakeETower(double vradial, 
                   double Tfreeze, 
                   double pTmin,
                   std::vector<particle> plist, 
                   std::vector<current> TypeOneSources,
                   std::vector<HadronizeCurrent> TypeTwoSources,
                   int coarse_level);
    void set_sigma(double _sigma){
        sigma = _sigma;
    }
    void FindJets(std::vector<double> Rs, 
                  double jetpTMin, 
                  double jetyMin, 
                  double jetyMax);
    void FindHF(std::vector<particle> plist);
    void CalcJetshape(std::vector<double> rbins);
    void LabelFlavor();
    void CorrHFET(std::vector<double> rbins);
    std::vector<Fjet> Jets;
    std::vector<particle> HFs;
    std::vector<Fjet> HFaxis;
private:
    int corp_index(double x, double xL, double xH, double dx, int Nx){
        if (x<xL || x>xH) return -1;
        else if (x==xH) return Nx-1;
        else return int((x-xL)/dx);
    };
    const int Neta, Nphi;
    const double etamax, etamin, phimax, phimin;
    const double deta, dphi;
    double sigma;
    MediumResponse MR;
    std::vector<std::vector<double> > PT;
    std::vector<std::vector<fourvec> > Pmu;
};

class LeadingParton{
   public:
   LeadingParton(std::vector<double> _pTbins);
   void add_event(std::vector<particle> plist, double sigma_gen);
   void write(std::string fheader);
   private:
   std::vector<double> pTbins, binwidth, nchg, npi, nD, nB;
   int NpT;
};

class JetStatistics{
   public:
   JetStatistics(std::vector<double> _pTbins, std::vector<double> Rs, std::vector<double> shapepTbins, std::vector<double> shaperbins);
   void add_event(std::vector<Fjet> jets, double sigma_gen);
   void write(std::string fheader);
   private:
   std::vector<double> pTbins, binwidth, shape_pTbins, shape_rbins, xJbins;
   std::vector<double> Rs, xJ;
   std::vector<std::vector<double> > shapes, Dshapes, Bshapes, dsigmadpT, dBdpT, dDdpT;
   int NpT, shape_NpT, shape_Nr;
};

class JetHFCorr{
   public:
   JetHFCorr(std::vector<double> pTHFbins, std::vector<double> rbins);
   void add_event(std::vector<Fjet> jets, std::vector<particle> HFs, 
                  double sigma_gen);
   void write(std::string fheader);
   private:
   std::vector<double> pTHFbins, rbins;
   std::vector<std::vector<double> > dDdr, dBdr;
};

class HFETCorr{
   public:
   HFETCorr(std::vector<double> pTHFbins, std::vector<double> rbins);
   void add_event(std::vector<Fjet> HFaxis, double sigma_gen);
   void write(std::string fheader);
   private:
   std::vector<double> pTHFbins, rbins;
   std::vector<std::vector<double> > D_dPTdr, B_dPTdr;
};


#endif
