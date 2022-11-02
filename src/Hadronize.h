#ifndef HADRONIZE_H
#define HADRONIZE_H

#include "Pythia8/Pythia.h"
#include "predefine.h"
#include <random>
using namespace Pythia8;

class HFdecayer{
public: 
    HFdecayer();
    int Decay(Pythia8::Particle pIn, std::vector<Pythia8::Particle> & pOut);
private:
    Pythia pythia;
};


class JetDenseMediumHadronize{
public:
   JetDenseMediumHadronize();
   int hadronize(std::vector<particle> partons, 
                 std::vector<particle> & hadrons, 
                 std::vector<particle> & thermal_partons,
                 double Q0, double Tf, int level);
private:
   Pythia pythia;
   HFdecayer Decay;
};



#endif
