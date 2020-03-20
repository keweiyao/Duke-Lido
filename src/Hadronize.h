// main21.cc is a part of the PYTHIA event generator.
// Copyright (C) 2018 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program.
// It illustrates how to feed in a single particle (including a resonance)
// or a toy parton-level configurations.

#include "Pythia8/Pythia.h"
#include "workflow.h"
#include <random>
using namespace Pythia8;

class JetDenseMediumHadronize{
public:
   JetDenseMediumHadronize();
   int hadronize(std::vector<particle> partons, 
                 std::vector<particle> & hadrons, 
                 std::vector<particle> & thermal_partons,
                 double Q0, int level);
private:
   Pythia pythia;
};
