#ifndef RANDOM_H
#define RANDOM_H
#include <random>
#include "lorentz.h"
namespace Srandom{
int getEnvSeed();
extern std::mt19937 gen;
extern std::uniform_real_distribution<double> sqrtZ;
extern std::uniform_real_distribution<double> rejection;
extern std::uniform_real_distribution<double> init_dis;
extern std::uniform_real_distribution<double> dist_phi;
extern std::uniform_real_distribution<double> dist_costheta;
extern std::normal_distribution<double> white_noise;
extern std::gamma_distribution<double> sample_E_over_T;
extern std::gamma_distribution<double> exp_dist;
int sample_flavor(int Nf);
bool binary_choice();
fourvec generate_thermal_parton_with_boost(double T, double vx, double vy, double vz);
}
#endif
