#include "Langevin.h"
#include "random.h"
#include "predefine.h"
#include "matrix_elements.h"
double const tiny = 1e-10;

double delta_qhat(int pid, double E, double M, double T){
	return 0.;
	/*
	int abspid = std::abs(pid);
        double EM = 1. + std::log (1+ (E-M)/qhat_params.b/T );
        double delta_qhat = qhat_params.K * std::pow(T, 3)
                    /(1. + std::pow(qhat_params.a*(T-Tc)/Tc, qhat_params.p))
                    /std::pow(EM, qhat_params.q);
            
        if (abspid !=21) {;
           return delta_qhat; 
        }
	else return delta_qhat*CA/CF;*/
}

double qhat_small_angle_LOpQCD(int pid, double E, double M, double T){
        double CR = (pid==21) ? CA : CF;
        double mD2 = t_channel_mD2->get_mD2(T);
        double Q2cut = cut*mD2;
        return alpha_s(Q2cut, T) * CR * T * mD2 * std::log(Q2cut/mD2);
}

double qhat_L_small_angle_LOpQCD(int pid, double E, double M, double T){
        double CR = (pid==21) ? CA : CF;
        double mD2 = t_channel_mD2->get_mD2(T);
        double Q2cut = cut*mD2;
        return alpha_s(Q2cut, T) * CR * T * .5*mD2 * std::log(Q2cut/mD2);
}


double qhat(int pid, double E, double M, double T){
	return  qhat_small_angle_LOpQCD(pid, E, M, T) 
	      + delta_qhat(pid, E, M, T);
}


double qhat_L(int pid, double E, double M, double T){
    return  qhat_L_small_angle_LOpQCD(pid, E, M, T);
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
	double p0 = std::sqrt(E0*E0 - M*M + 1e-6);

	double kt = qhat(pid, E0, M, T)/2.;
        double kl = qhat_L(pid, E0, M, T);
	double minf = std::sqrt(t_channel_mD2->get_mD2(T)/2.);
	double Ed = std::max(E0, minf);
	double drag = kl/(2.*Ed*T);
        double damp = 1.-drag*dt;		   
	double Ct = std::sqrt(kt*dt);
	double Cl = std::sqrt(kl*dt);
   
        pOut.a[1] = Ct * Srandom::white_noise(Srandom::gen);
        pOut.a[2] = Ct * Srandom::white_noise(Srandom::gen);
        pOut.a[3] = p0 * damp 
                  + Cl * Srandom::white_noise(Srandom::gen);
        pOut.a[0] = std::sqrt(M*M + std::pow(pOut.x(),2) 
		  	+ std::pow(pOut.y(),2) + std::pow(pOut.z(),2) );

	// rotate back to the original frame
	pOut = pOut.rotate_back(pIn_cell);
	// boost back to lab frame
	pOut = pOut.boost_back(v[0], v[1], v[2]);
}

