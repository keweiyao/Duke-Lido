#include "Xsection.h"

template<size_t N>
Xsection<N>::Xsection(std::string Name, 
					boost::property_tree::ptree config):
StochasticBase<N>(Name, config)
{
}

template<size_t N>
void Xsection<N>::init(void){
	std::cout << "init table" << std::endl;
}

template<size_t N>
void Xsection<N>::compute(void){
	std::cout << "computing table" << std::endl;
}

template<size_t N>
void Xsection<N>::sample(std::vector<double> arg, 
						std::vector< std::vector<double> > & FS){
	std::cout << "sampling table" << std::endl;
}

template class Xsection<2>;
