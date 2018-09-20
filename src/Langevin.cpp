#include "Langevin.h"
#include "random.h"
#include "predefine.h"
#include "matrix_elements.h"
double A=0., B=0.;
double const tiny = 1e-10;
double logof2 = std::log(2);
// for quarks, upto t=mD^2

double qhat_pQCD(int pid, double E, double T){
	double factor = 1.0;
	if (pid==21) factor = CA/CF; 
	double alphas = alpha_s(0, T);
	double mD2 = t_channel_mD2->get_mD2(T);
	return A*factor*alphas*CF*T * mD2;
}
double dqhat_pQCD_dp2(int pid, double E, double T){
	return 0.;
}

/*
double kperp(double E, double M, double T){
	return std::pow(T,3)*( A + B/(E*T) );
}

double kpara(double E, double M, double T){
	return std::pow(T,3)*( A + B/(E*T) );
}

double dkpara_dp2(double E, double M, double T){
	return std::pow(T,3)*B*(-1.)/(2.*E*E*E*T);
}*/


void initialize_transport_coeff(double _A, double _B){
	A = _A; B = _B;
	std::cout << "A = " << A << ", B = " << B << std::endl;
};

void Ito_update(int pid, double dt_lab, double M, double T, std::vector<double> v, 
						const fourvec & pIn, fourvec & pOut){
	// Boost pIn to medium frame
	auto pIn_cell = pIn.boost_to(v[0], v[1], v[2]);
	// Boost dt to medium frame
	double dt = dt_lab*pIn_cell.t()/pIn.t();
	// imaging rotating to a frame where pIn lies on z-axis
	double E0 = pIn_cell.t();
	double p0 = std::sqrt(E0*E0 - M*M+1e-9);

	double kt = qhat_pQCD(pid, E0, T)/2.;
	double kl = kt;
	double dkl_dp2 = dqhat_pQCD_dp2(pid, E0, T)/2.;
	double drag = kl/(2.*E0*T) - (kl - kt)/std::pow(p0, 2) - dkl_dp2;
		   
	double Ct = std::sqrt(kt*dt);
	double Cl = std::sqrt(kl*dt);

    pOut.a[1] = Ct * Srandom::white_noise(Srandom::gen);
    pOut.a[2] = Ct * Srandom::white_noise(Srandom::gen);
    pOut.a[3] = p0 * (1. - drag * dt) + Cl * Srandom::white_noise(Srandom::gen);
    pOut.a[0] = std::sqrt(M*M + std::pow(pOut.x(),2) 
					+ std::pow(pOut.y(),2) + std::pow(pOut.z(),2) );

	// rotate back to the original frame
	pOut = pOut.rotate_back(pIn_cell);
	// boost back to lab frame
	pOut = pOut.boost_back(v[0], v[1], v[2]);
}

