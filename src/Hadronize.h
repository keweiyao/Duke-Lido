#ifndef HADRONIZE_H
#define HADRONIZE_H

#include "Pythia8/Pythia.h"
#include "predefine.h"
#include <random>
using namespace Pythia8;

class JetDenseMediumHadronize{
public:
   JetDenseMediumHadronize();
   int hadronize(std::vector<particle> partons, 
                 std::vector<particle> & hadrons, 
                 std::vector<particle> & thermal_partons,
                 double Q0, double Tf, int level);
private:
   Pythia pythia;
};

class HFHadronize{
public:
   HFHadronize();
   int hadronize(particle pIn, std::vector<particle> & pOut, double Q0, double Tf);
private:
   Pythia pythia;
};


#endif
