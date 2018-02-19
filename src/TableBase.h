#ifndef TABLE_BASE_H
#define TABLE_BASE_H

#include <vector>
#include <string>
#include <boost/multi_array.hpp>
#include <iostream>
#include "lorentz.h"

typedef std::vector<double> Dvec;
typedef std::vector<size_t> Svec;

// Base class of a table of type T with dimension N
template <typename T, size_t N>
class TableBase{
protected:
    const std::string _Name;
    const size_t _rank, _power_rank;
    Svec _shape;
    Dvec _low, _high;
    Dvec _step;
    boost::multi_array<T, N> _table;
public:
	TableBase(std::string, Svec, Dvec, Dvec);
	T InterpolateTable(Dvec values);
    void SetTableValue(Svec index, T v);
    bool Save(std::string);
    bool Load(std::string);
};





#endif

