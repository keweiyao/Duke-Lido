#ifndef Xsection_H
#define Xsection_H

#include <cstdlib>
#include <vector>
#include <string>
#include <random>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include "StochasticBase.h"


template <const char * str, size_t N, typename F>
class Xsection: public virtual StochasticBase<N> {
private:
    scalar find_max(std::vector<double> parameters);
	scalar calculate_scalar(std::vector<double> parameters);
	fourvec calculate_fourvec(std::vector<double> parameters);
	tensor calculate_tensor(std::vector<double> parameters);
	std::vector<double> _IS_masses;
	std::vector<double> _FS_masses;
	std::vector<int> _IS_types;
	std::vector<int> _FS_types;
        int _process_id;
	F _f;// the matrix element
public:
	Xsection(std::string Name, std::string configfile, F f);
	bool sample(std::vector<double> arg,
                        int incoming_hard_pid,
			std::vector<fourvec> & FS,
                        std::vector<int> & pids );
};

#endif
