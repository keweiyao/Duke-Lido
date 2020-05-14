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
    std::vector<double> shape, Leadingshape, dndr;
    std::vector<particle> Ftags;
    double sigma;
};
void TestSource(
             MediumResponse MR,
             std::vector<current> jlist,
             std::vector<HadronizeCurrent> slist,
             std::string fname
     );
std::vector<Fjet> FindJetTower(MediumResponse MR,
             std::vector<particle> plist,
             std::vector<current> SourceList,
             std::vector<HadronizeCurrent> HadronizeList,
             std::vector<double> Rs,
	     std::vector<double> rbins,
             double jetpTMin,
             double jetyMin,
             double jetyMax,
             double sigma_gen,
	     double pTmin);
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
   std::vector<Fjet> AllJets;
   std::vector<double> pTbins, binwidth, shape_pTbins, shape_rbins, xJbins;
   std::vector<double> shape_w, Rs, xJ;
   std::vector<std::vector<double> > shapes, dnchdr, dsigmadpT, dBdpT, dDdpT, dB0dpT, dD0dpT;
   int NpT, shape_NpT, shape_Nr;
};

#endif
