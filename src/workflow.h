#ifndef WORKFLOW_H
#define WORKFLOW_H

#include <boost/variant/variant.hpp>
#include "Rate.h"
#include "Onium_Rate.h"
#include <vector>
#include <map>

struct particle{
	// mass, x, p, t, all in units of [GeV^a]
	int pid;
	double mass, weight;
	bool is_vac, is_virtual;
	
	double T0, mfp0, Tf; // production temperature, local mfp
	fourvec x0; // production location
	fourvec x; // current location
	fourvec p0; // production momentum
	fourvec p; // current momentum

	std::vector<particle> radlist;
	std::vector<double> vcell;

	void freestream(double dt){
		double a = dt/p.t();
		x.a[0] = x.t() + dt;
		x.a[1] = x.x() + p.x()*a;
		x.a[2] = x.y() + p.y()*a;
		x.a[3] = x.z() + p.z()*a;
	}
};

typedef Rate<LO, 2, 2, double(*)(const double, void*)> Rate22;
typedef Rate<GB, 2, 2, double(*)(const double*, void*)> Rate23;
typedef Rate<GB, 2, 4, double(*)(const double*, void*)> Rate32;
typedef EffRate12<2, double(*)(const double*, void*)> Rate12;
typedef EffRate21<2, double(*)(const double*, void*)> Rate21;
typedef OniumDissoRate22<2, double(*)(double, void *)> OniumD22;
typedef OniumRecoRate22<3, double(*)(double *, std::size_t, void*)> OniumR22;
typedef OniumDissoRate23q<2, 
  double(*)(double *, std::size_t, void *), 
  double(*)(double, void *), 
  double(*)(double, double)> OniumD23q;
typedef OniumRecoRate32q<3, 
  double(*)(double *, std::size_t, void *), 
  double(*)(double, void *), 
  double(*)(double, double)> OniumR32q;
typedef boost::variant<Rate22,   // case 0
                       Rate23,   // case 1
                       Rate32,   // case 2
                       Rate12,   // case 3 
                       Rate21,   // case 4
                       OniumD22, // case 5 
                       OniumR22,  // case 6
                       OniumD23q, // case 7
                       OniumR32q // case 8
               > Process;
extern std::map<int, std::vector<Process>> AllProcesses;

void init_process(Process& r, std::string mode, std::string table_path);
void initialize(std::string mode, std::string setting_path, std::string table_path, double mu, double afix, 
		double K, double a, double b, double p, double q, double gamma, double cut, double Rvac);

int OneBodyUpdate_Parton(double dt, double temp, std::vector<double> v3cell, particle & pIn, std::vector<particle> & pOut_list);
int OneBodyUpdate_Onium(double dt, double temp, std::vector<double> v3cell, particle & pIn, std::vector<particle> & pOut_list);
int TwoBodyUpdate_QQbar(double dt_lab, double temp, std::vector<double> v3cell, 
		particle & P1, particle & P2, std::vector<particle> & pOut_list);

double formation_time(fourvec p, fourvec k, double T, int split);
double calcualte_dt_from_dtau(fourvec x, fourvec p, double tau, double dtau);
void output(const std::vector<particle> plist, std::string fname);
void output_oscar(const std::vector<particle> plist, int abspid, std::string fname);
double mean_pT(const std::vector<particle> plist);
double mean_E(const std::vector<particle> plist);
#endif
