#ifndef JET_FINDING_H
#define JET_FINDING_H

#include <fstream>
#include <vector>
#include <string>
#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/FastJet3.h"
#include "lorentz.h"
#include "workflow.h"

struct current{
    fourvec p;
    double chetas, shetas;
    double cs;
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

void TestSource(
             std::vector<current> jlist,
             std::string fname);
void FindJetTower(std::vector<particle> plist, 
             std::vector<current> jlist,
             std::vector<double> Rs,
             double jetpTMin, 
             double jetyMin, 
             double jetyMax,
             std::string fname, 
             double sigma_gen);

#endif
