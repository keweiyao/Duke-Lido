#ifndef LIDO_H
#define LIDO_H

#include "predefine.h"
#include "collision_manager.h"
#include "Hadronize.h"
#include <vector>


// for 1->2, 2->2 and 2->3 collisions
class lido{
private:
    std::unique_ptr<collision_manager> CM;
    double Lido_Ecut;
    int FrameChoice;
    void Diffusion(particle & p, double dt, double T, std::vector<double> v);
public:
    lido(std::string setting_path, std::string table_path, 
         std::vector<double> parameters);
    int update_single_particle(
        double dt, double T, std::vector<double> v3,
        particle & pIn, std::vector<particle> & pOut_list);
    void FreeStream(particle & p, double dt);
    void set_frame(int _choice){
        FrameChoice = _choice;
    }; // Cartisian / Bjorken
};


#endif
