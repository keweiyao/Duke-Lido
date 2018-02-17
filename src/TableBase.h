#ifndef TABLE_BASE_H
#define TABLE_BASE_H

#include <cstdlib>
#include <vector>
#include <string>
#include <random>
#include <stdarg.h>
#include <boost/multi_array.hpp>
#include <iostream>
#include <type_traits>
#include "utility.h"


// Base class of a table of type T with dimension N
template <typename T, size_t N>
class TableBase{
protected:
    const std::string _Name;
    const size_t _rank, _power_rank;
    const std::vector<size_t> _shape;
    boost::multi_array<T, N> _table;
    T (*_approximating_function)(std::vector<double> inputs);
public:
	TableBase(std::string, T (*f)(std::vector<double>), std::vector<size_t>);
	T InterpolateTable(std::vector<double> values);
    void SetTableValue(std::vector<size_t> index, T v);
    void echo(void);
};





#endif

