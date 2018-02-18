#include "TableBase.h"
#include <iostream>

template <typename T, size_t N>
TableBase<T, N>::TableBase(std::string Name, T (*f)(Dvec), Svec shape, Dvec low, Dvec high):
_Name(Name), _rank(N), _power_rank(std::pow(2, _rank)), 
_shape(shape), _low(low), _high(high),
_table(_shape), _approximating_function(f) {
	std::cout<<_Name << " dim=" << _rank <<std::endl;
	for(auto i=0; i<_rank; ++i){
		_step.push_back((high[i]-low[i])/(shape[i]-1));
	}
	for(auto& item : _step)	std::cout << item << " "; 
	std::cout<<std::endl;
}

template <typename T, size_t N>
T TableBase<T, N>::InterpolateTable(Dvec values){
   Svec start_index;
   Dvec w;
   for(auto i=0; i<_rank; ++i) {
       auto x = (values[i]-_low[i])/_step[i];
       size_t nx = size_t(std::floor(x));
       double rx = x-nx;
       start_index.push_back(nx);
       w.push_back(rx);
   }
   Svec index(_rank);
   T result{0.};
   for(auto i=0; i<_power_rank; ++i) {
        auto W = 1.0;
        for (auto j=0; j<_rank; ++j) {
            index[j] = start_index[j] + ((i & ( 1 << j )) >> j);
            W *= (index[j]==start_index[j])?(1.-w[j]):w[j];
        }
        result = result + _table(index)*W;
   }
   return result;
}

template <typename T, size_t N>
void TableBase<T, N>::SetTableValue(Svec index, T v){
    _table(index) = v;
}


template class TableBase<double, 2>;
template class TableBase<double, 3>;
template class TableBase<double, 4>;
template class TableBase<fourvec, 2>;
template class TableBase<fourvec, 3>;


