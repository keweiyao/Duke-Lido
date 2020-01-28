#include "random.h"
#include <string>
#include <iostream>

namespace Srandom{
int getEnvSeed(){
    char *val = std::getenv("lseed");
    int seed = (val == NULL) ? -1 : std::stoi(std::string(val));
    std::cout << "Lido seed = " << seed;
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
std::normal_distribution<double> white_noise(0.0, 1.0);
int sample_flavor(int Nf){
    int flavor = (gen()%Nf) + 1;
    if (gen()%2 == 0) 
        return flavor;
    else 
        return -flavor;
}

}
