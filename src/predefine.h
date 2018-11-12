#ifndef PREDEFINE_H
#define PREDEFINE_H

#include <cmath>
#include <H5Cpp.h>
#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>

extern char LO[];
extern char GB[];

extern bool type1_warned;
extern bool type2_warned;
extern bool type3_warned;

//=============useful constants=============================================
extern const double c4d9;
extern const double c1d9;
extern const double c16pi;
extern const double c48pi;
extern const double c16pi2;
extern const double c64d9pi2;
extern const double c256pi4;
extern const double fmc_to_GeV_m1;
// number of color=3 (3*3-1 = 8 gluons), number of flavor=3, (u,d,s quark)
extern const int Nc, nf;
extern const double CF;
extern const double CA;
extern const double CF_over_CA;

// the prefractor for gluon debye mass with Boltzmann statistics
// mD^2 = 8\pi*(Nc+nf)*alpha_s*T^2 ~ 15*alpha_s*T^2
// Note that if using quantum statistics, this will be
// mD^2 = 4\pi/3*(Nc+nf/2)*alpha_s*T^2 ~ 18*alpha_s*T^2
extern const double pf_g; // prefractor for gluon self energy^2

// For QCD coupling constant
// alpha_s = alpha_0 = 4\pi/(11Nc/3-2Nf/3)/log(Q^2/LambdaQCD^2)
extern const double alpha0; // alpha_s(Q2 = e*Lambda2)
extern const double Lambda; // [GeV] Lambda QCD = 0.2 GeV
extern const double Lambda2; // [GeV^2] Lambda QCD squared
extern const double mu2_left; // minimum cut on Q2, where alpha = alpha_0

// helper function for read/write hdf5 scalar attributes
template <typename T> inline const H5::PredType& type();
template <> inline const H5::PredType& type<size_t>() { return H5::PredType::NATIVE_HSIZE; }
template <> inline const H5::PredType& type<double>() { return H5::PredType::NATIVE_DOUBLE; }
template <> inline const H5::PredType& type<int>() { return H5::PredType::NATIVE_INT; }

template <typename T>
void hdf5_add_scalar_attr(
  const H5::Group& gp, const std::string& name, const T& value) {
  const auto& datatype = type<T>();
  auto attr = gp.createAttribute(name.c_str(), datatype, H5::DataSpace{});
  attr.write(datatype, &value);
}

template <typename T>
void hdf5_read_scalar_attr(
  const H5::Group& gp, const std::string& name, T& value) {
  const auto& datatype = type<T>();
  auto attr = gp.openAttribute(name.c_str());
  attr.read(datatype, &value);
}

// Fortran style float format
class ff {
public:
    ff(double x): value(x) {}
    const double value;
	friend std::ostream & operator<< (std::ostream & stream, const ff & x) {
		// So that the log does not scream
		if (x.value == 0.) {
		    stream << "0.000000D+00";
		    return stream;
		}
		int exponent = floor(log10(std::abs(x.value)));
		double base = x.value / pow(10, exponent);
		// Transform here
		base /= 10;
		exponent += 1;
		if (base >= 0){
			std::stringstream buff;
			buff << std::setw(8) << std::setprecision(6);
			buff << std::fixed << base;
			stream << std::setw(8) << buff.str();
		}
		else{
			std::stringstream buff;
			buff << std::setw(8) << std::setprecision(6);
			buff << std::fixed << base;
			//Checking if we have a leading minus sign
			std::string newbase = "-" + buff.str().substr(2, buff.str().size()-1);
			stream << std::setw(8) << newbase;
		}
		
		if (exponent >= 0) stream << "D+" << std::setw(2) << std::setfill('0') << exponent;
		else stream << "D-" << std::setw(2) << std::setfill('0') << std::abs(exponent);
		return stream;
	}
};

#endif
