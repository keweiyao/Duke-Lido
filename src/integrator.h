#ifndef INTEGRATOR_H
#define INTEGRATOR_H
#include <iostream>
#include <cmath>

#include <functional>
#include <memory>
#include <utility>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

/* Modified from here 
@MISC {27248,
    TITLE = {C++ library for numerical intergration (quadrature)},
    AUTHOR = {Henri Menke (https://scicomp.stackexchange.com/users/24680/henri-menke)},
    HOWPUBLISHED = {Computational Science Stack Exchange},
    NOTE = {URL:https://scicomp.stackexchange.com/q/27248 (version: 2017-06-27)},
    EPRINT = {https://scicomp.stackexchange.com/q/27248},
    URL = {https://scicomp.stackexchange.com/q/27248}
}
A wrapper aournd the terrible GSL interface...
*/

template < typename F >
class gsl_quad_1d{
  F f;
  int limit;
  std::unique_ptr < gsl_integration_workspace,
                    std::function < void(gsl_integration_workspace*) >
                    > workspace;

  static double gsl_wrapper(double x, void * p)
  {
    gsl_quad_1d * t = reinterpret_cast<gsl_quad_1d*>(p);
    return t->f(x);
  }

public:
  gsl_quad_1d(F f, int limit): 
  f(f), limit(limit), 
  workspace(gsl_integration_workspace_alloc(limit), gsl_integration_workspace_free)
  {}

  double integrate(double min, double max, double epsabs, double epsrel, double &error)
  {
    gsl_function gsl_f;
    gsl_f.function = &gsl_wrapper;
    gsl_f.params = this;

    double result;
    if ( !std::isinf(min) && !std::isinf(max) ){
      gsl_integration_qags ( &gsl_f, min, max,
                             epsabs, epsrel, limit,
                             workspace.get(), &result, &error );
    }
    else if ( std::isinf(min) && !std::isinf(max) ){
      gsl_integration_qagil( &gsl_f, max,
                             epsabs, epsrel, limit,
                             workspace.get(), &result, &error );
    }
    else if ( !std::isinf(min) && std::isinf(max) ){
      gsl_integration_qagiu( &gsl_f, min,
                             epsabs, epsrel, limit,
                             workspace.get(), &result, &error );
    }
    else{
      gsl_integration_qagi ( &gsl_f,
                             epsabs, epsrel, limit,
                             workspace.get(), &result, &error );
    }
    return result;
  }
};

template < typename F >
double quad_1d(F func,
            std::pair<double,double> const& range, double &error,
            double epsabs = 1.e-4, double epsrel = 1.e-4,
            int limit = 1000){
  return gsl_quad_1d<F>(func, limit).integrate(range.first, range.second, epsabs, epsrel, error);
}

#endif
