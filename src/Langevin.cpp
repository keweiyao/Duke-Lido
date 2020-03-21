#include "Langevin.h"
#include "random.h"
#include "predefine.h"
#include "matrix_elements.h"
double const tiny = 1e-10;

double delta_qhat(int pid, double E, double M, double T){
        double CR = (pid==21) ? CA : CF;
        int absid = std::abs(pid);
        double EM = 1. + std::log (1+ (E-M)/qhat_params.b/T );
        double delta_qhat = CR/CF*qhat_params.K * std::pow(T, 3)
                /(1. + std::pow(qhat_params.a*(T-Tc)/Tc, qhat_params.p))
                /std::pow(EM, qhat_params.q);
        return delta_qhat;
}

double qhat_small_angle_LOpQCD(int pid, double E, double M, double T){
        double CR = (pid==21) ? CA : CF;
        double mD2 = t_channel_mD2->get_mD2(T);
        double Q2cut = std::min(cut*mD2, 6*E*T);
        return alpha_s(Q2cut, T) * CR * T * mD2 * std::log(1.+Q2cut/mD2);
}

double qhat_L_small_angle_LOpQCD(int pid, double E, double M, double T){
        double CR = (pid==21) ? CA : CF;
        double minf2 = .5*t_channel_mD2->get_mD2(T);
        double Q2cut = std::min(cut*minf2, 3*E*T);
        return alpha_s(Q2cut, T) * CR * T * minf2 * std::log(1.+Q2cut/minf2);
}


double qhat(int pid, double E, double M, double T){
	return  qhat_small_angle_LOpQCD(pid, E, M, T) 
	      + delta_qhat(pid, E, M, T);
}


double qhat_L(int pid, double E, double M, double T){
	double m0;
	double minf = std::sqrt(t_channel_mD2->get_mD2(T)/2.);
    m0 = std::max(minf, M);
    return  qhat_L_small_angle_LOpQCD(pid, E, M, T)
              + delta_qhat(pid, E, m0, T)/2. * std::pow(E/m0, qhat_params.gamma);                       
}

double dqhat_L_dp2(int pid, double E, double M, double T){
	double p2 = E*E - M*M + tiny;
	double dp2 = p2*.05;
	double Eprime = std::sqrt(E*E + dp2);
	return (qhat_L(pid, Eprime, M, T) - qhat_L(pid, E, M, T) ) /dp2;
}

void Ito_update(int pid, double dt_lab, double M, double T, std::vector<double> v, 
						const fourvec & pIn, fourvec & pOut){
	// Boost pIn to medium frame
	auto pIn_cell = pIn.boost_to(v[0], v[1], v[2]);
	// Boost dt to medium frame
	double dt = dt_lab*pIn_cell.t()/pIn.t();
	// imaging rotating to a frame where pIn lies on z-axis
	double E0 = pIn_cell.t();
	double p0 = std::sqrt(E0*E0 - M*M+1e-9);

	double kt = qhat(pid, E0, M, T)/2.;
        double kl = qhat_L(pid, E0, M, T);
	double dkl_dp2 = dqhat_L_dp2(pid, E0, M, T);
	double drag = kl/(2.*E0*T);
                   // - (kl - kt)/std::pow(p0, 2) - dkl_dp2;
		   
	double Ct = std::sqrt(kt*dt);
	double Cl = std::sqrt(kl*dt);

        pOut.a[1] = Ct * Srandom::white_noise(Srandom::gen);
        pOut.a[2] = Ct * Srandom::white_noise(Srandom::gen);
        pOut.a[3] = p0 * (1. - drag * dt) 
                  + Cl * Srandom::white_noise(Srandom::gen);
        pOut.a[0] = std::sqrt(M*M + std::pow(pOut.x(),2) 
		  	+ std::pow(pOut.y(),2) + std::pow(pOut.z(),2) );

	// rotate back to the original frame
	pOut = pOut.rotate_back(pIn_cell);
	// boost back to lab frame
	pOut = pOut.boost_back(v[0], v[1], v[2]);
}

