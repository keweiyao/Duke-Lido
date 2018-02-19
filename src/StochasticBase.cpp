#include "StochasticBase.h"

template<size_t N>
StochasticBase<N>::StochasticBase(std::string Name)
: _Name(Name)
{
	std::cout << __func__ << " " << _Name << std::endl;
}

