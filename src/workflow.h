#ifndef WORKFLOW_H
#define WORKFLOW_H

#include <boost/variant/variant.hpp>
#include "Rate.h"
#include <vector>
#include <map>

// particle data type
struct pregluon{
	fourvec p0, k1, kn;
	double t0, T0;
	double local_mfp;
	bool is_vac;
};

struct particle{
	// mass, x, p, t, all in units of [GeV^a]
	int pid;
	bool freezeout;
	double mass;
	fourvec x;
	fourvec p;
	std::vector<pregluon> radlist, abslist;
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
	double weight;
};

typedef Rate<LO, 2, 2, double(*)(const double, void*)> Rate22;
typedef Rate<GB, 2, 2, double(*)(const double*, void*)> Rate23;
typedef Rate<GB, 2, 4, double(*)(const double*, void*)> Rate32;
typedef EffRate12<2, double(*)(const double*, void*)> Rate12;
typedef EffRate21<2, double(*)(const double*, void*)> Rate21;
typedef boost::variant<Rate22, Rate23, Rate32, Rate12, Rate21> Process;
extern std::map<int, std::vector<Process>> AllProcesses;

void init_process(Process& r, std::string mode, std::string table_path);
void initialize(std::string, std::string setting_path, std::string table_path, double mu, double alpha_s_fixed, double A, double B);

int gluon_elastic_scattering(double dt, double temp, std::vector<double> v3cell, fourvec incomping_p, fourvec & outgoing_p);
int update_particle_momentum_Lido(double dt, double temp, std::vector<double> v3cell, particle & pIn);

double formation_time(fourvec p, fourvec k, double M, double T);
double calcualte_dt_from_dtau(fourvec x, fourvec p, double tau, double dtau);
void output_oscar(const std::vector<particle> plist, std::string fname);
double mean_pT(const std::vector<particle> plist);

#endif
