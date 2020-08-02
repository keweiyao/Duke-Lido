#ifndef APPROX_FUNCTIONS_H
#define APPROX_FUNCTIONS_H
#include <vector>
#include "lorentz.h"

// Xsection
scalar approx_X22(std::vector<double> params);
scalar approx_X22QQbar(std::vector<double> params);
scalar approx_dX22_max(std::vector<double> params);
#endif
