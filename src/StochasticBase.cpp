#include <cstdlib>
#include <vector>
#include <string>
#include <random>
#include <boost/multi_array.hpp>
#include "StochasticBase.h"
//	double (*dXdPS)(double * PS, size_t n_dims, void * params);
// this base defines how a stochasitc object
// 1) calculate the probablity of a event with input parameters
// 2) tabulate with I/O of the probablity
// 3) interpolate the probablity over input parameters
// 4) given inputs, sample the output parameters

class StochasticBase{
protected:
    std::random_device _rd;
    std::mt19937 _gen;
	virtual void _tabulate_prob() = 0;
	virtual void _save_to_file() = 0;
	virtual void _read_from_file() = 0;
    std::string _Name;
public:
	StochasticBase(std::string Name);
	virtual double interpolate_prob(double * arg) = 0; 
	virtual double calculate_prob(double * arg) = 0;
	virtual void sample_output(double * arg, std::vector< std::vector<double> > & FS) = 0;
};


StochasticBase::StochasticBase(std::string Name)
: Name(Name)
{
	std::cout << __func__ << " " << _Name << std::endl;
}

