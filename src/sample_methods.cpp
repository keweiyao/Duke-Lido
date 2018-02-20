#include "sample_methods.h"
#include <random>
#include <fstream>
#include <algorithm>

void rejection_1d::build_interval(double xL, double xH, double fxL, double fxH){
	double xM = 0.5*(xL+xH);
	double fxM = f(&xM, 1, params);
	double mid_h = 0.5*(fxL+fxH);
	if (fxM < mid_h*0.1){
		build_interval(xL, xM, fxL, fxM);
		build_interval(xM, xH, fxM, fxH);
	}
	else{
		rectangle A; A.xL = xL; A.dx = xH-xL; A.fL = fxL; A.df = fxH-fxL; A.w = A.dx*mid_h;
		intervals.push_back(A);
		total_weight += A.w;
		return;
	}
}

double rejection_1d::sample(double (*f_) (double * x, size_t n_dims, void * params), double xlo_, double xhi_, void * params_){
	f = f_;
	xlo = xlo_;
	xhi = xhi_;
	params = params_;
	intervals.clear();
	total_weight = 0.0;
			
	double fL = f(&xlo, 1, params), fH = f(&xhi, 1, params);
	build_interval(xlo, xhi, fL, fH);
	double cumulate = 0.0;
	for (auto&& ele : intervals){
		cumulate += ele.w/total_weight;
		ele.w = cumulate;
	}
	double r1, r2, xl=0., dx=0., fl=0., df=0., fh=0., lambda, xtry;
	do{
		r1 = std::rand()*1./RAND_MAX;
		for (auto&& ele : intervals){
			if (ele.w > r1) {
				xl = ele.xL; dx = ele.dx;
				fl = ele.fL; df = ele.df;
				fh = fl+df;
				break;
			}
		}
		r2 = std::rand()*1./RAND_MAX;
		lambda = ( std::sqrt(fl*fl*(1.-r2) + fh*fh*r2) -fl )/df;
		xtry = xl + lambda*dx;
		}while (f(&xtry, 1, params)/(fl + df*lambda)*RAND_MAX < std::rand());
	return xtry;
}

double rejection_1d::plain_sample(double (*f_) (double * x, size_t n_dims, void * params), double xlo_, double xhi_, void * params_){
	f = f_;
	xlo = xlo_;
	xhi = xhi_;
	params = params_;
	double dx = xhi-xlo, xtry;
	double fmax = std::max(f(&xlo, 1, params), f(&xhi, 1, params));
	do{
		xtry = xlo + std::rand()*dx/RAND_MAX;
	}while (f(&xtry, 1, params)/fmax*RAND_MAX < std::rand());
	return xtry;
}

// ----------Affine-invariant metropolis sample-------------------
AiMS::AiMS(void)
:	a(0.5), rd(), gen(rd()), sqrtZ(std::sqrt(1./a), std::sqrt(a)),
	reject(0.0, 1.0)
{
}
void AiMS::initialize(void){
    std::uniform_real_distribution<double> init_dis(0, 1);
	for (size_t i=0; i<Nwalker; ++i){
		do{
			for (size_t j=0; j < n_dims; ++j) walkers[i].posi[j] = guessl[j] + (guessh[j]-guessl[j])*init_dis(gen);
			walkers[i].P = f(walkers[i].posi, n_dims, params);
		} while(walkers[i].P <= 1e-22);
		for (size_t j=0; j < n_dims; ++j)
				buff_walkers[i].posi[j] = walkers[i].posi[j];
		buff_walkers[i].P = walkers[i].P;
	}
}
void AiMS::update(void){
	size_t ri;
	double sqz, z, Ptry, Paccept;
	double * xtry = new double[n_dims];
	walker w, wr;
    for (size_t i=0; i<Nwalker; ++i){
		do{ 
			ri = rd() % Nwalker;
		}while(i==ri);
		w = walkers[i];
		wr = walkers[ri];
		sqz = sqrtZ(gen);
		z = sqz*sqz;
		for (size_t j=0; j < n_dims; ++j) xtry[j] = wr.posi[j] + z*(w.posi[j] - wr.posi[j]);
		Ptry = f(xtry, n_dims, params);
		Paccept = Ptry/w.P*std::pow(z, n_dims-1);
		if (Paccept >= 1.0){
			for (size_t j=0; j < n_dims; ++j) buff_walkers[i].posi[j] = xtry[j];
			buff_walkers[i].P = Ptry;
		}
		else if (Paccept >= reject(gen)){
			for (size_t j=0; j < n_dims; ++j) buff_walkers[i].posi[j] = xtry[j];
			buff_walkers[i].P = Ptry;
		}
	}
	for (size_t i=0; i<Nwalker; ++i){
		for (size_t j=0; j < n_dims; ++j) walkers[i].posi[j] = buff_walkers[i].posi[j];
		walkers[i].P = buff_walkers[i].P;
	}
	delete[] xtry;
}

std::vector<double> AiMS::sample(double (*f_) (double*, size_t, void*), size_t n_dims_, void * params_, double * guessl_, double * guessh_){
	walkers.clear();
	f = f_; n_dims = n_dims_; params = params_; guessl = guessl_; guessh = guessh_;
	Nwalker = n_dims*4;
	walkers.resize(Nwalker);
	buff_walkers.resize(Nwalker);
	for (auto&& w : walkers) w.posi = new double[n_dims];
	for (auto&& w : buff_walkers) w.posi = new double[n_dims];

	initialize();
	for (size_t i = 0; i<Nwalker*20; i++) { update(); }

	std::vector<double> result;
	result.resize(n_dims);
	for (size_t i = 0; i<n_dims; i++) result[i] = walkers[0].posi[i];

	for (auto&& w : walkers) delete[] w.posi;
	for (auto&& w : buff_walkers) delete[] w.posi;

	return result;
}

