#ifndef StochasticBase_H
#define StochasticBase_H

#include <cstdlib>
#include <vector>
#include <string>
#include <random>
#include <boost/multi_array.hpp>
#include "sample_methods.h"
//	double (*dXdPS)(double * PS, size_t n_dims, void * params);
// this base defines how a stochasitc object
// 1) calculate the probablity of a event with input parameters (0th moment)
// 2) calculate the first and second moments for the output parameters
// 3) tabulate with I/O of the probablity
// 4) interpolate the probablity over input parameters
// 5) given inputs, sample the output parameters

template <size_t input_dim>
class StochasticBase{
protected:
    std::random_device rd;
    std::mt19937 gen;
	virtual void tabulate_prob() = 0;
	virtual void save_to_file() = 0;
	virtual void read_from_file() = 0;
public:
	Xsection(double (*dXdPS_)(double *, size_t, void *), double M1_, std::string name_, bool refresh);
	// arg = [s, T] fot X22, arg = [s, T, dt] for X23, arg = [s, T, s1k, s2k] for f32
	virtual double interpolate_prob(double * arg) = 0; 
	virtual double calculate_prob(double * arg) = 0;
	virtual void sample_output(double * arg, std::vector< std::vector<double> > & FS) = 0;
};

#endif
