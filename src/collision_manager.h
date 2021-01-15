#ifndef Collision_Manager_H
#define Collision_Manager_H

#include <boost/variant/variant.hpp>
#include "Rate.h"
#include <vector>
#include <map>

typedef Rate<HS2PP, 2, 2, double(*)(const double, void*)> Rate22;
typedef Rate<HS2PPP, 2, 2, double(*)(const double*, void*)> Rate23;
typedef EffRate12<2, double(*)(const double*, void*)> Rate12;
typedef boost::variant<Rate22, Rate23, Rate12> Process;

// for 1->2, 2->2 and 2->3 collisions
class collision_manager{
private:
    std::map<int, std::vector<Process> > AllProcesses;
    void init_process(Process &r, std::string mode, std::string table_path);
    double get_effective_mfp_soft(
                              double Ecell, double x,
                              int idA, int idB, int idC,
                              double Temp
                              );
    double get_effective_mfp_hard(
                              double Ecell, double x,
                              int idA, int idB, int idC,
                              double Temp
                              );
public:
    collision_manager(std::string mode, std::string setting_path,
                    std::string table_path, std::vector<double> parameters);
    int sample(particle & pIn, 
              std::vector<particle> & pOut_list,
              double Temp, 
              std::vector<double> v3cell,
              double dt_lab);
    // parameters = mu, afix, theta, cut
};


#endif
