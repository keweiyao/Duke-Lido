#ifndef TABLE_BASE_H
#define TABLE_BASE_H

#include <cstdlib>
#include <vector>
#include <string>
#include <random>
#include <stdarg.h>
#include <boost/multi_array.hpp>
#include <iostream>

struct ScalarData{
    double m;
    friend inline std::ostream& operator<<(std::ostream& os, const ScalarData& a) {
      os << a.m;
      return os;
    }
};

struct VectorData{
    double t, x, y, z;
    friend inline std::ostream& operator<<(std::ostream& os, const VectorData& a) {
      os << "x^\\mu = " << a.t << ' ' << a.x << ' '
                                  << a.y << ' ' << a.z;
      return os;
    }
};

struct TensorData{
    double tt, tx, ty, tz;
    double xt, xx, xy, xz;
    double yt, yx, yy, yz;
    double zt, zx, zy, zz;
    friend inline std::ostream& operator<<(std::ostream& os, const TensorData& a) {
      os << "x^\\mu\\nu = " 
         << a.tt << ' ' << a.tx << ' ' << a.ty << ' ' << a.tz << std::endl 
         << a.xt << ' ' << a.xx << ' ' << a.xy << ' ' << a.xz << std::endl 
		 << a.yt << ' ' << a.yx << ' ' << a.yy << ' ' << a.yz << std::endl 
		 << a.zt << ' ' << a.zx << ' ' << a.zy << ' ' << a.zz;
      return os;
    }
};

template <typename T, size_t N>
class TableBase{
protected:
    std::string _Name;
    const size_t _dims;
    T (*_approximating_function)(std::vector<double> inputs);
    boost::multi_array<T, N> _table;
public:
	TableBase(std::string Name, T (*f)(std::vector<double>));
	virtual T InterpolateTable(size_t, ...) = 0;
    virtual void SetTableAt(size_t, ...) = 0;
    virtual T GetTableAt(size_t, ...) = 0;
};

// 2D table
class ScalarTable2D: public TableBase<ScalarData, 2> {
public:
	ScalarTable2D(std::string, ScalarData (*f)(std::vector<double>));
	ScalarData InterpolateTable(size_t, ...);
	void SetTableAt(size_t, ...);
	ScalarData GetTableAt(size_t, ...);
};

class VectorTable2D: public TableBase<VectorData, 2> {
public:
	VectorTable2D(std::string, VectorData (*f)(std::vector<double>));
	VectorData InterpolateTable(size_t, ...);
	void SetTableAt(size_t, ...);
	VectorData GetTableAt(size_t, ...);
};

class TensorTable2D: public TableBase<TensorData, 2> {
public:
	TensorTable2D(std::string, TensorData (*f)(std::vector<double>));
	TensorData InterpolateTable(size_t, ...);
	void SetTableAt(size_t, ...);
	TensorData GetTableAt(size_t, ...);
};

// 3D table
class ScalarTable3D: public TableBase<ScalarData, 3> {
public:
	ScalarTable3D(std::string, ScalarData (*f)(std::vector<double>));
	ScalarData InterpolateTable(size_t, ...);
	void SetTableAt(size_t, ...);
	ScalarData GetTableAt(size_t, ...);
};

// 4D table
class ScalarTable4D: public TableBase<ScalarData, 4> {
public:
	ScalarTable4D(std::string, ScalarData (*f)(std::vector<double>));
	ScalarData InterpolateTable(size_t, ...);
	void SetTableAt(size_t, ...);
	ScalarData GetTableAt(size_t, ...);
};

template <size_t N> ttype();

#endif

