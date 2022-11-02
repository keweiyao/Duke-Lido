#ifndef JET_FINDING_H
#define JET_FINDING_H

#include <fstream>
#include <vector>
#include <string>
#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/FastJet3.h"
#include "lorentz.h"
#include "TableBase.h"
#include "predefine.h"

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
    MyInfo(int ieta, int iphi) : _ieta(ieta), _iphi(iphi) {}; 
    int ieta() const {return _ieta;}
    int iphi() const {return _iphi;}
    private:
    int _ieta, _iphi;
};
struct Fjet{
    fourvec pmu;
    double R, pT, M2, phi, eta;
    std::vector<double> shape, dNdz, dDdz, dBdz, dNdpT;
    int flavor;
    double sigma;
};
void TestSource(
             MediumResponse MR,
             std::vector<current> jlist,
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
                   int coarse_level,
		   bool charged_jet);
    void set_sigma(double _sigma){
        sigma = _sigma;
    }
    void FindJets(std::vector<double> Rs, 
                  double jetpTMin, 
                  double jetyMin, 
                  double jetyMax, 
		  bool charged_trigger);
    void FindHF(std::vector<particle> plist);
    void CalcJetshape(std::vector<double> rbins);
    void Frag(std::vector<double> zbins, std::vector<double>  zpTbins);
    void LabelFlavor();
    std::vector<Fjet> Jets;
    std::vector<particle> HFs;
    std::vector<std::vector<double> > PT;
private:
    int corp_index(double x, double xL, double xH, double dx, int Nx){
        if (x<xL || x>xH) return -1;
        else if (x==xH) return Nx-1;
        else return int((x-xL)/dx);
    };
    const int Neta, Nphi;
    const double etamax, etamin, phimax, phimin;
    const double deta, dphi;
    double sigma, vradial, Tfreeze;
    std::vector<particle> plist;
    std::vector<current> clist;
    std::vector<double> Rs;
    MediumResponse MR;
    std::vector<std::vector<fourvec> > Pmu;
    Pythia8::Pythia HF2mu;
    std::map<int,double> BRs;
};

class LeadingParton{
   public:
   LeadingParton();
   void add_event(std::vector<particle> plist, double sigma_gen);
   void write(std::string fheader);
   private:
   std::vector<double> pTbins, binwidth, nchg, npi, nstrange, nD, nB;
   std::vector<std::vector<double> > v2chg, v2strange, v2pi, v2D, v2B,
                        v3chg, v3strange, v3pi, v3D, v3B;
   int NpT;
};

class JetStatistics{
   public:
   JetStatistics(std::vector<double> Rs, 
      std::vector<double> shaperbins, 
      std::vector<double> Fragszbins, 
      std::vector<double> FragszpTbins);
   void add_event(std::vector<Fjet> jets, double sigma_gen, fourvec x0);
   void write(std::string fheader);
   private:
   std::vector<double> pTbins, binwidth, 
	   shape_pTbins, shape_rbins, 
	   xJbins, xJ_pTbins, xJ_W, 
	   Frag_pTbins, Frag_zbins, Frag_zpTbins, 
           Frags_W, Frags_D_W, Frags_B_W,
           Shape_W, Shape_D_W, Shape_B_W,
           leading_yield, subleading_yield;
   std::vector<std::vector<std::vector<double> > > JQ2, JQ3, JQ4,  D_in_jet, B_in_jet;
   std::vector<double> Rs;
   std::vector<std::vector<double> > shapes, xJ, 
	   Dshapes, Bshapes, 
	   dsigmadpT, dBdpT, dDdpT,
	   Frags, Frags_pT, 
           Frags_D, Frags_D_pT,
           Frags_B, Frags_B_pT,
           DijetInfo,  D_in_jet_W, B_in_jet_W;
   int NpT, shape_NpT, shape_Nr;
};

class JetHFCorr{
   public:
   JetHFCorr(std::vector<double> rbins);
   void add_event(std::vector<Fjet> jets, std::vector<particle> HFs, 
                  double sigma_gen);
   void write(std::string fheader);
   private:
   std::vector<double> pTHFbins, rbins, dDdr_W, dBdr_W;
   std::vector<std::vector<double> > dDdr, dBdr;
};

#endif
