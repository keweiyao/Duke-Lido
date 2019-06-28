#ifndef Onium_Rate_H
#define Onium_Rate_H

#include <cstdlib>
#include <vector>
#include <string>
#include <random>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include "StochasticBase.h"
#include "predefine.h"
#include "Onium_predefine.h"

// Quarkonium 2->2 disso rate: QQbar[nl] + g --> Q + Qbar
// N: dimension of rate table: 2: velocity(v), tempature(T)
// F: diff rate function type
// What it does: given v, T, determine rate, and sample final state
template <size_t N, typename F>
class OniumDissoRate22: public virtual StochasticBase<N>{
private:
    scalar find_max(std::vector<double> parameters);
	scalar calculate_scalar(std::vector<double> parameters);
	fourvec calculate_fourvec(std::vector<double> parameters);
	tensor calculate_tensor(std::vector<double> parameters);
	F _f; // the kernel
	double _mass, _Enl, _aB; // by convention _Enl>0, mass is 2*M*Mbar/(M+Mbar)
    int _n, _l;
	bool _active;
public:
	OniumDissoRate22(std::string Name, std::string configfile, int n, int l, F f);
	void sample(std::vector<double> arg, 
				std::vector< fourvec > & FS);
	bool IsActive(void) {return _active;}
    int n(){return _n;}
    int l(){return _l;}
    double aB(){return _aB;}
    double Enl(){return _Enl;}
    void FlavorContent(int & Q, int & Qbar){
       if (_mass<1.8) {
           Q = 4;
           Qbar = -4;
       }
       else if (_mass<4) {
           Q = 5;
           Qbar = -4;
       }
       else{
           Q = 5;
           Qbar = -5;
       }
    }
};

// Quarkonium 2->2 reco rate: Q + Qbar --> QQbar[nl] + g
// N: dimension of rate table: 3: velocity(v), tempature(T). relative momentum (p)
// F: diff rate function type
template <size_t N, typename F>
class OniumRecoRate22: public virtual StochasticBase<N>{
private:
    scalar find_max(std::vector<double> parameters);
	scalar calculate_scalar(std::vector<double> parameters);
	fourvec calculate_fourvec(std::vector<double> parameters);
	tensor calculate_tensor(std::vector<double> parameters);
	F _f; // the kernel
	double _mass, _Enl, _aB; // by convention _Enl>0
    int _n, _l;
	bool _active;
public:
	OniumRecoRate22(std::string Name, std::string configfile, int n, int l, F f);
	void sample(std::vector<double> arg, 
				std::vector< fourvec > & FS);
	bool IsActive(void) {return _active;}
    int n(){return _n;}
    int l(){return _l;}
    double aB(){return _aB;}
    double Enl(){return _Enl;}
};


// Quarkonium 2->3 quark disso rate: QQbar[nl] + q --> Q + Qbar + q
// N: dimension of rate table: 2: velocity(v), tempature(T)
// F: diff rate function type
template <size_t N, typename F>
class OniumDissoRate23q: public virtual StochasticBase<N>{
private:
    scalar find_max(std::vector<double> parameters);
    scalar calculate_scalar(std::vector<double> parameters);
    fourvec calculate_fourvec(std::vector<double> parameters);
    tensor calculate_tensor(std::vector<double> parameters);
    F _f1; // the kernel
    F _f2; // the kernel
    double _mass, _Enl, _aB; // by convention _Enl>0, mass is 2*M*Mbar/(M+Mbar)
    int _n, _l;
    bool _active;
public:
    OniumDissoRate23q(std::string Name, std::string configfile, int n, int l, F f1, F f2);
    void sample(std::vector<double> arg,
                std::vector< fourvec > & FS);
    bool IsActive(void) {return _active;}
    int n(){return _n;}
    int l(){return _l;}
    void FlavorContent(int & Q, int & Qbar){
        if (_mass<1.8) {
            Q = 4;
            Qbar = -4;
        }
        else if (_mass<4) {
            Q = 5;
            Qbar = -4;
        }
        else{
            Q = 5;
            Qbar = -5;
        }
    }
};

// Quarkonium 3->2 quark reco rate: Q + Qbar + q --> QQbar[nl] + q
// N: dimension of rate table: 3: velocity(v), tempature(T), relative momentum(p)
// F: diff rate function type
template <size_t N, typename F>
class OniumRecoRate32q: public virtual StochasticBase<N>{
private:
    scalar find_max(std::vector<double> parameters);
    scalar calculate_scalar(std::vector<double> parameters);
    fourvec calculate_fourvec(std::vector<double> parameters);
    tensor calculate_tensor(std::vector<double> parameters);
    F _f1; // the kernel
    F _f2; // the kernel
    double _mass, _Enl, _aB; // by convention _Enl>0, mass is 2*M*Mbar/(M+Mbar)
    int _n, _l;
    bool _active;
public:
    OniumRecoRate32q(std::string Name, std::string configfile, int n, int l, F f1, F f2);
    void sample(std::vector<double> arg,
                std::vector< fourvec > & FS);
    bool IsActive(void) {return _active;}
    int n(){return _n;}
    int l(){return _l;}
};

#endif
