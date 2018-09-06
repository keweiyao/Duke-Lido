#ifndef WORKFLOW_H
#define WORKFLOW_H

#include <boost/variant/variant.hpp>
#include "Rate.h"
#include <vector>
#include <map>

struct pregluon{
	fourvec p0, k1, kn;
	double t0, T0;
	double local_mfp;
};

struct particle{
	// mass, x, p, t, all in units of [GeV^a]
	int pid;
	bool freezeout;
	double mass;
	fourvec x;
	fourvec p;
	std::vector<pregluon> radlist, abslist;
	double t_rad, t_abs;
	fourvec p0;
	std::vector<double> vcell;
	double Tf;
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
typedef boost::variant<Rate22, Rate23, Rate32, Rate12, Rate21> Process;

extern std::map<int, std::vector<Process>> AllProcesses;
void initialize(std::string, std::string path, double mu, double alpha_s_fixed, double A, double B);

double formation_time(fourvec p, fourvec k, double M, double T, double mpf);

int gluon_elastic_scattering(double dt, double temp, std::vector<double> v3cell, fourvec incomping_p, fourvec & outgoing_p);

//typedef Rate<GBHT, 3, 3, double(*)(const double*, void*)> Rate23;
//typedef Rate<GBHT, 3, 4, double(*)(const double*, void*)> Rate32;
int update_particle_momentum_HT(double dt, double temp, std::vector<double> v3cell, particle & pIn);

int update_particle_momentum_Lido(double dt, double temp, std::vector<double> v3cell, particle & pIn);

std::vector<double> probe_test(double M, double E0, double T, double dt, int Nsteps,
				int Nparticles, std::string mode, double mu, double const_alphas, double A, double B);

std::vector<double> Bjorken_test(double M, double E0, double T0, double t0, double dt, int Nsteps, int Nparticles, std::string mode, double mu, double const_alphas, double A, double B);

#endif
