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

typedef std::vector<double> Dvec;
typedef std::vector<size_t> Svec;

// Base class of a table of type T with dimension N
template <typename T, size_t N>
class TableBase{
protected:
    const std::string _Name;
    const size_t _rank, _power_rank;
    const Svec _shape;
    const Dvec _low, _high, _step;
    boost::multi_array<T, N> _table;
    T (*_approximating_function)(Dvec inputs);
public:
	TableBase(std::string, T (*f)(Dvec), Svec);
	T InterpolateTable(Dvec values);
    void SetTableValue(Svec index, T v);
    void echo(void);
};





#endif

