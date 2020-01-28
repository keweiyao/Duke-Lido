#include <fstream>
#include <string>
#include <math.h>
#include <vector>
#include "jet_finding.h"
#include "workflow.h"
#include "lorentz.h"

void FindJet(std::vector<particle> plist, 
             std::vector<fourvec> jlist,
             double jetRadius, 
             double jetpTMin, 
             double jetyMin, 
             double jetyMax,
             std::string fname, double sigma_gen) {
    int power = -1; // power=-1: anti-kT
    //fastjet.init();
    // Define a jet
    fastjet::JetDefinition jetDef(fastjet::genkt_algorithm, jetRadius, power);
    // Objetcts to operate on
    std::vector<fastjet::PseudoJet> fjInputs;
    // cuts
    fastjet::Selector select_rapidity = fastjet::SelectorRapRange(jetyMin, jetyMax);
    for (auto & p : plist){
        fastjet::PseudoJet fp(p.p.x(), p.p.y(), p.p.z(), p.p.t());
        fp.set_user_info(new MyInfo(p.pid, p.origin, sigma_gen));
        fjInputs.push_back(fp);
    }
    for (auto & p : jlist){
        fastjet::PseudoJet fp(p.x(), p.y(), p.z(), p.t());
        fp.set_user_info(new MyInfo(0, 3, sigma_gen));
        fjInputs.push_back(fp);
    }
    
    std::vector<fastjet::PseudoJet> inclusiveJets, jets;
    fastjet::ClusterSequence clustSeq(fjInputs, jetDef);
    inclusiveJets = clustSeq.inclusive_jets(jetpTMin);
    std::vector<fastjet::PseudoJet> selected_jets = select_rapidity(inclusiveJets);
    jets = sorted_by_pt(selected_jets);
    std::ofstream f(fname.c_str(), std::ios_base::app);
    for (auto & j : jets){
        double pT = std::sqrt(j.px()*j.px() + j.py()*j.py());
        double pabs = std::sqrt(pT*pT + j.pz()*j.pz());
        double phi = std::atan2(j.py(), j.px());
        double M2 = j.e()*j.e() - pabs*pabs;
        double rap = .5*std::log((pabs+j.pz())/(pabs-j.pz()));
        f << pT << " " << phi << " " << rap << " " << M2 << " " 
          << sigma_gen << " "  
          << j.e() << " "  << j.px() << " " << j.py() << " " << j.pz() 
          << std::endl;
    }
}

