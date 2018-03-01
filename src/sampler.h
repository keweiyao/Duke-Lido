#ifndef SAMPLER_H
#define SAMPLER_H
#include <iostream>
#include <cmath>

#include <functional>
#include <memory>
#include <utility>
#include <random>


template < typename F >
double sample_1d(F f, std::pair<double,double> const& range, double fmax){
  	double y, x, xlow=range.first, xhigh=range.second;
  	double interval = xhigh-xlow;
	do{
		x = xlow+std::rand()*interval/RAND_MAX;
		y = f(x)/fmax;
	}while(std::rand()*1./RAND_MAX>y);
	return x;
}

template < typename F >
std::vector<double> sample_nd(F f, int dim, std::vector<std::pair<double,double>> const& range, double fmax){
	double * x = new double[dim];
  	double y;
  	double * interval = new double[dim];
	int counter = 0;
	for(int i=0; i<dim; i++) interval[i] = range[i].second - range[i].first;
	do{
		// random choice
		for(int i=0; i<dim; i++) 
			x[i] = range[i].first+std::rand()*interval[i]/RAND_MAX;
		y = f(x)/fmax;
		counter ++;
	}while(std::rand()*1./RAND_MAX>y && counter < 10000);
	std::vector<double> res(dim);
	for(int i=0; i<dim; i++) res[i] = x[i];
	delete[] x;
	delete[] interval;
	std::cout << "try = " << counter<< std::endl;
	return res;
}


#endif 
