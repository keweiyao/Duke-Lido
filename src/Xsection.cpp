#include "Xsection.h"


template<size_t N>
Xsection<N>::Xsection(std::string Name, 
					boost::property_tree::ptree config):
StochasticBase<N>(Name, config)
{
}

template<size_t N>
void Xsection<N>::sample(std::vector<double> arg, 
						std::vector< std::vector<double> > & FS){
	std::cout << "sampling table" << std::endl;
}

template<size_t N>
scalar Xsection<N>::calculate_scalar(std::vector<double> parameters){
	scalar res;
	res.s = parameters[0]*parameters[1];
	return res;
}

template class Xsection<2>;
template class Xsection<3>;
template class Xsection<4>;
