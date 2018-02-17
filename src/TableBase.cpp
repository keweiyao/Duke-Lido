#include "TableBase.h"
#include <iostream>

Dvec calculate_steps(Dvec low, Dvec high, Dvec)

template <typename T, size_t N>
TableBase<T, N>::TableBase(std::string Name, T (*f)(Dvec), Svec shape, Dvec low, Dvec high):
_Name(Name), _rank(N), _power_rank(std::pow(2, _rank)), 
_shape(shape), _low(low), _high(high), _step(calculate_steps)
_table(_shape), _approximating_function(f) {
	std::cout<<_Name << " dim=" << _rank <<std::endl;
}

template <typename T, size_t N>
T TableBase<T, N>::InterpolateTable(Dvec values){
   Svec index(_rank);
   Dvec w = values;
   T result{0};
   for(auto i=0; i<_power_rank; ++i) {
        auto W = 1.0;
        for (auto j=0; j<_rank; ++j) {
            index[j] = (i & ( 1 << j )) >> j;
            W *= (index[j]==0)?(1-w[j]):w[j];
        }
        result = result + _table(index)*W;
        std::cout << result << ' ';
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
template class TableBase<VectorT, 2>;
template class TableBase<TensorT, 2>;


