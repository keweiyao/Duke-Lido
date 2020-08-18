#include "random.h"
#include <string>
#include <iostream>
#include "lorentz.h"
namespace Srandom{
int getEnvSeed(){
    char *val = std::getenv("lseed");
    int seed = (val == NULL) ? -1 : std::stoi(std::string(val));
    std::cout << "Lido seed = " << seed << std::endl;
    return seed;
}
std::mt19937 gen(  (getEnvSeed() < 0) ? std::random_device{}() : getEnvSeed());
const double AMC = 4.0;
std::uniform_real_distribution<double> sqrtZ(
       std::sqrt(1./AMC), std::sqrt(AMC)
    );
std::uniform_real_distribution<double> rejection(0.0, 1.0);
std::uniform_real_distribution<double> init_dis(0.0, 1.0);
std::uniform_real_distribution<double> dist_phi(0.0, 2.0*M_PI);
std::uniform_real_distribution<double> dist_costheta(-1.0, 1.0);
std::gamma_distribution<double> sample_E_over_T(3.0, 1.0);
std::gamma_distribution<double> exp_dist(1.0,1.0);
std::normal_distribution<double> white_noise(0.0, 1.0);
int sample_flavor(int Nf){
    int flavor = (gen()%Nf) + 1;
    if (gen()%2 == 0) 
        return flavor;
    else 
        return -flavor;
}
bool binary_choice(){
    if (gen()%2==0) return true;
    else return false;
}
fourvec generate_thermal_parton_with_boost(double T, double vx, double vy, double vz){
    // randomly sample the four momentum
    double se;
    do {
        se = T*Srandom::sample_E_over_T(Srandom::gen);
    } while (se>5.*T);
    double phi = Srandom::dist_phi(Srandom::gen);
    double cos = Srandom::dist_costheta(Srandom::gen);
    double sin = std::sqrt(1.-cos*cos);
    fourvec p{se, se*sin*std::cos(phi), se*sin*std::sin(phi),se*cos};
    return p.boost_back(vx, vy, vz);
}
}
