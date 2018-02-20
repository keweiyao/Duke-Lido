#ifndef SAMPLE_METHODS_H
#define SAMPLE_METHODS_H

#include <vector>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <random>

struct rectangle{
	double xL, dx;
	double fL, df;
	double w;
};

class rejection_1d{
private:
	double (*f) (double * x, size_t n_dims, void * params);
	void * params;
	double xlo, xhi, total_weight;
	std::vector<rectangle> intervals;
	void build_interval(double xL, double xH, double fxL, double fxH);
public:
	rejection_1d(void){};
	double sample(double (*f_) (double * x, size_t n_dims, void * params), double xlo_, double xhi_, void * params_);
	double plain_sample(double (*f_) (double * x, size_t n_dims, void * params), double xlo_, double xhi_, void * params_);
};


// ----------Affine-invariant metropolis sample-------------------
struct walker{
	double * posi;
	double P;
};

class AiMS{
private:
	double (*f) (double*, size_t, void*);	
    void * params;
	size_t n_dims, Nwalker;
	std::vector<walker> walkers, buff_walkers;
	void initialize(void);
	void update(void);
	double a;
	std::random_device rd;
    std::mt19937 gen;
    std::uniform_real_distribution<double> sqrtZ;
	std::uniform_real_distribution<double> reject;
	double * guessl, * guessh;
public:
	AiMS(void);
	std::vector<double> sample(double (*f_) (double*, size_t, void*), size_t n_dims_, void * params_, double * guessl_, double * guessh_);
};

#endif
