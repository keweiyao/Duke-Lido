#include "workflow.h"
#include "simpleLogger.h"
#include <fstream>
#include "random.h"
#include "matrix_elements.h"
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/variant/get.hpp>
#include <boost/any.hpp>
#include <boost/foreach.hpp>
#include "logo.h"
#include "Langevin.h"
#include "predefine.h"

std::map<int, std::vector<Process>> AllProcesses;
std::ofstream ofs ("hq_gluon_accept.dat", std::ofstream::app);

void init_process(Process &r, std::string mode, std::string table_path)
{
	switch (r.which())
	{
	case 0:
		if (boost::get<Rate22>(r).IsActive())
			if (mode == "new")
			{
				boost::get<Rate22>(r).initX(table_path);
				boost::get<Rate22>(r).init(table_path);
			}
			else
			{
				boost::get<Rate22>(r).loadX(table_path);
				boost::get<Rate22>(r).load(table_path);
			}
		else
			return;
		break;
	case 1:
		if (boost::get<Rate23>(r).IsActive())
			if (mode == "new")
			{
				boost::get<Rate23>(r).initX(table_path);
				boost::get<Rate23>(r).init(table_path);
			}
			else
			{
				boost::get<Rate23>(r).loadX(table_path);
				boost::get<Rate23>(r).load(table_path);
			}
		else
			return;
		break;
	case 2:
		if (boost::get<Rate32>(r).IsActive())
			if (mode == "new")
			{
				boost::get<Rate32>(r).initX(table_path);
				boost::get<Rate32>(r).init(table_path);
			}
			else
			{
				boost::get<Rate32>(r).loadX(table_path);
				boost::get<Rate32>(r).load(table_path);
			}
		else
			return;
		break;
	case 3:
		if (boost::get<Rate12>(r).IsActive())
			if (mode == "new")
			{
				boost::get<Rate12>(r).init(table_path);
			}
			else
			{
				boost::get<Rate12>(r).load(table_path);
			}
		else
			return;
		break;
	case 4:
		if (boost::get<Rate21>(r).IsActive())
			if (mode == "new")
			{
				boost::get<Rate21>(r).init(table_path);
			}
			else
			{
				boost::get<Rate21>(r).load(table_path);
			}
		else
			return;
		break;
	default:
		exit(-1);
		break;
	}
}

void initialize(std::string mode, std::string setting_path, std::string table_path)
{
	print_logo();

	boost::property_tree::ptree config;
	std::ifstream input(setting_path);
	read_xml(input, config);

	double mu=config.get("mu", 2.0); //alphas(max(Q, mu*pi*T)), active when afix < 0
	double afix=config.get("afix", -1.0);  // fixed coupling, negative value for running coupling
	double K=config.get("K", 0.0);
	double a=config.get("a", 1.0);
	double b=config.get("b", 1.0);
	double p=config.get("p", 1.0);
	double q=config.get("q", 1.0);
	double gamma=config.get("gamma", 1.0); // K, a, b, p, q, gamma are used in non-pert. qhat parametrization
										   // Its effect is off when K = 0.
	double cut=config.get("cut", 1.0); // Separation scale between diffusion
										 // and scattering both at leading order (weak coupled)
	double Rvac=config.get("Rvac", 1.0); // Vaccum veto region parameter

	Lido_Ecut = config.get("Ecut", 2.0);

	initialize_mD_and_scale(0, mu, afix, cut, Rvac);
	initialize_transport_coeff(K, a, b, p, q, gamma);
	echo();

	// for gluon
	AllProcesses[21] = std::vector<Process>();
	AllProcesses[21].push_back(Rate22("Boltzmann/gq2gq", setting_path, dX_gq2gq_dt));	  // 2->2, index = 0
	AllProcesses[21].push_back(Rate22("Boltzmann/gg2gg", setting_path, dX_gg2gg_dt));	  // 2->2, index = 1
	AllProcesses[21].push_back(Rate23("Boltzmann/gq2gqg", setting_path, M2_gq2gqg));	   // 2->3, index = 2
	AllProcesses[21].push_back(Rate23("Boltzmann/gg2ggg", setting_path, M2_gg2ggg));	   // 2->3, index = 3
	AllProcesses[21].push_back(Rate32("Boltzmann/gqg2gq", setting_path, M2_gqg2gq));	   // 3->2, index = 4
	AllProcesses[21].push_back(Rate32("Boltzmann/ggg2gg", setting_path, M2_ggg2gg));	   // 3->2, index = 5
	AllProcesses[21].push_back(Rate12("Boltzmann/g2gg", setting_path, LGV_g2gg));		   // 1->2, index = 6
	AllProcesses[21].push_back(Rate21("Boltzmann/gg2g", setting_path, LGV_gg2g));		   // 2->1, index = 7
	AllProcesses[21].push_back(Rate23("Boltzmann/gq2qqqbar", setting_path, M2_gq2qqqbar)); // 2->3, index = 8
	AllProcesses[21].push_back(Rate23("Boltzmann/gg2qgqbar", setting_path, M2_gg2qgqbar)); // 2->3, index = 9
	AllProcesses[21].push_back(Rate12("Boltzmann/g2qqbar", setting_path, LGV_g2qqbar));	// 1->2, index = 10

	// for u, d, s flavor
	AllProcesses[123] = std::vector<Process>();
	AllProcesses[123].push_back(Rate22("Boltzmann/qq2qq", setting_path, dX_Qq2Qq_dt)); // 2->2, index = 0
	AllProcesses[123].push_back(Rate22("Boltzmann/qg2qg", setting_path, dX_Qg2Qg_dt)); // 2->2, index = 1
	AllProcesses[123].push_back(Rate23("Boltzmann/qq2qqg", setting_path, M2_Qq2Qqg));  // 2->3, index = 2
	AllProcesses[123].push_back(Rate23("Boltzmann/qg2qgg", setting_path, M2_Qg2Qgg));  // 2->3, index = 3
	AllProcesses[123].push_back(Rate32("Boltzmann/qqg2qq", setting_path, M2_Qqg2Qq));  // 3->2, index = 4
	AllProcesses[123].push_back(Rate32("Boltzmann/qgg2qg", setting_path, M2_Qgg2Qg));  // 3->2, index = 5
	AllProcesses[123].push_back(Rate12("Boltzmann/q2qg", setting_path, LGV_q2qg));	 // 1->2, index = 6
	AllProcesses[123].push_back(Rate21("Boltzmann/qg2q", setting_path, LGV_qg2q));	 // 2->1, index = 7

	// for charm flavor
	AllProcesses[4] = std::vector<Process>();
	AllProcesses[4].push_back(Rate22("Boltzmann/cq2cq", setting_path, dX_Qq2Qq_dt)); // 2->2, index = 0
	AllProcesses[4].push_back(Rate22("Boltzmann/cg2cg", setting_path, dX_Qg2Qg_dt)); // 2->2, index = 1
	AllProcesses[4].push_back(Rate23("Boltzmann/cq2cqg", setting_path, M2_Qq2Qqg));  // 2->3, index = 2
	AllProcesses[4].push_back(Rate23("Boltzmann/cg2cgg", setting_path, M2_Qg2Qgg));  // 2->3, index = 3
	AllProcesses[4].push_back(Rate32("Boltzmann/cqg2cq", setting_path, M2_Qqg2Qq));  // 3->2, index = 4
	AllProcesses[4].push_back(Rate32("Boltzmann/cgg2cg", setting_path, M2_Qgg2Qg));  // 3->2, index = 5
	AllProcesses[4].push_back(Rate12("Boltzmann/c2cg", setting_path, LGV_q2qg));	 // 1->2, index = 6
	AllProcesses[4].push_back(Rate21("Boltzmann/cg2c", setting_path, LGV_qg2q));	 // 2->1, index = 7

	// for bottom flavor
	AllProcesses[5] = std::vector<Process>();
	AllProcesses[5].push_back(Rate22("Boltzmann/bq2bq", setting_path, dX_Qq2Qq_dt)); // 2->2, index = 0
	AllProcesses[5].push_back(Rate22("Boltzmann/bg2bg", setting_path, dX_Qg2Qg_dt)); // 2->2, index = 1
	AllProcesses[5].push_back(Rate23("Boltzmann/bq2bqg", setting_path, M2_Qq2Qqg));  // 2->3, index = 2
	AllProcesses[5].push_back(Rate23("Boltzmann/bg2bgg", setting_path, M2_Qg2Qgg));  // 2->3, index = 3
	AllProcesses[5].push_back(Rate32("Boltzmann/bqg2bq", setting_path, M2_Qqg2Qq));  // 3->2, index = 4
	AllProcesses[5].push_back(Rate32("Boltzmann/bgg2bg", setting_path, M2_Qgg2Qg));  // 3->2, index = 5
	AllProcesses[5].push_back(Rate12("Boltzmann/b2bg", setting_path, LGV_q2qg));	 // 1->2, index = 6
	AllProcesses[5].push_back(Rate21("Boltzmann/bg2b", setting_path, LGV_qg2q));	 // 2->1, index = 7

	// Initialzie all processes
	BOOST_FOREACH (Process &r, AllProcesses[123])
		init_process(r, mode, table_path);
	BOOST_FOREACH (Process &r, AllProcesses[4])
		init_process(r, mode, table_path);
	BOOST_FOREACH (Process &r, AllProcesses[5])
		init_process(r, mode, table_path);
	BOOST_FOREACH (Process &r, AllProcesses[21])
		init_process(r, mode, table_path);
}

// split=1: q->q+g, colors = 1 - x + CF/CA * x^2
// split=2: g->g+g, colors = 1 - x + x^2
// split=3: g->q+qbar, colors = 1 - x*CA/CF + x^2*CA/CF
double formation_time(fourvec p, fourvec k, double T, int split)
{
	double E0 = p.t();
	double x = k.t() / E0;
	double mg2 = t_channel_mD2->get_mD2(T) / 2.,
		   mass_sqrs = 0.0, colors = 1.;
	if (split == 1)
	{
		// q --> q + g
		colors = 1. - x + CF / CA * x * x;
		mass_sqrs = x * x * dot(p, p) + (1. - x) * mg2;
	}
	if (split == 2)
	{
		// g --> g + g
		colors = 1. - x + x * x;
		mass_sqrs = (1. - x + x * x) * mg2;
	}
	if (split == 3)
	{
		// g --> q + qbar
		colors = 1. - x * CA / CF + x * x * CA / CF;
		mass_sqrs = (1 - CA / CF * x + CA / CF * x * x) * mg2;
	}

	double QT2 = measure_perp(p, k).pabs2() * colors;
	double tauf = 2 * x * (1 - x) * E0 / (QT2 + mass_sqrs);
	return tauf;
}

int update_particle_momentum_Lido(
	double dt, double temp, std::vector<double> v3cell,
	particle &pIn, std::vector<particle> &pOut_list)
{  
	auto p00 = pIn.p;
	//const double Lido_Ecut=5.*std::max(temp, 0.154);
	//const double Lido_Ecut = 2.0;


	pOut_list.clear();
	pIn.Tf = temp;
	pIn.vcell[0] = v3cell[0]; 
	pIn.vcell[0] = v3cell[1]; 
	pIn.vcell[0] = v3cell[3];
	auto x0 = pIn.x;

	double tau0 = std::sqrt(std::pow(pIn.x.t(),2) - std::pow(pIn.x.z(),2));
	pIn.freestream(dt);
	double tau1 = std::sqrt(std::pow(pIn.x.t(),2) - std::pow(pIn.x.z(),2));
    double dtau = tau1 - tau0;

	//ignore soft partons with a constant momentum cut
	if (pIn.p.boost_to(v3cell[0], v3cell[1], v3cell[2]).t() < Lido_Ecut || temp < 0.154)
	{
		pIn.radlist.clear();
		pOut_list.push_back(pIn);
		return pOut_list.size();
	}

	//LOG_INFO << "pid: " << pIn.pid;
	//LOG_INFO << pIn.p;
	//LOG_INFO << pIn.x;

	//transform light quark pid to 123
	//transform back at the end
	int temp_pid = pIn.pid;
	if (std::abs(pIn.pid) <= 3) pIn.pid = 123;


	int absid = std::abs(pIn.pid);

	// Apply diffusion and update particle momentum
	fourvec pnew;
	Ito_update(pIn.pid, dt, pIn.mass, temp, v3cell, pIn.p, pnew);
	pIn.p = pnew;

	// Apply large angle scattering, and diffusion induced radiation
	auto p_cell = pIn.p.boost_to(v3cell[0], v3cell[1], v3cell[2]);
	double dt_cell = dt / pIn.p.t() * p_cell.t();
	double E_cell = p_cell.t();

	// For diffusion induced radiation, qhat_Soft is an input
	double qhatg = qhat(21, E_cell, 0., temp);

	// Total rate for the
	int channel = 0;
	std::vector<double> P_channels(AllProcesses[absid].size());
	double P_total = 0., dR;

	BOOST_FOREACH (Process &r, AllProcesses[absid])
	{
		switch (r.which())
		{
		case 0:
			if (boost::get<Rate22>(r).IsActive())
				dR = boost::get<Rate22>(r).GetZeroM({E_cell, temp}).s;
			else
				dR = 0.0;
			P_channels[channel] = P_total + dR * dt_cell;
			break;
		case 1:
			if (boost::get<Rate23>(r).IsActive() && (!pIn.is_virtual))
				dR = boost::get<Rate23>(r).GetZeroM({E_cell, temp}).s;
			else
				dR = 0.0;
			P_channels[channel] = P_total + dR * dt_cell;
			break;
		case 2:
			if (boost::get<Rate32>(r).IsActive() && (!pIn.is_virtual))
				dR = boost::get<Rate32>(r).GetZeroM({E_cell, temp}).s;
			else
				dR = 0.0;
			P_channels[channel] = P_total + dR * dt_cell;
			break;
		case 3:
			if (boost::get<Rate12>(r).IsActive() && (!pIn.is_virtual))
				dR = qhatg * boost::get<Rate12>(r).GetZeroM({E_cell, temp}).s;
			else
				dR = 0.0;
			P_channels[channel] = P_total + dR * dt_cell;
			break;
		case 4:
			if (boost::get<Rate21>(r).IsActive() && (!pIn.is_virtual))
				dR = qhatg * boost::get<Rate21>(r).GetZeroM({E_cell, temp}).s;
			else
				dR = 0.0;
			P_channels[channel] = P_total + dR * dt_cell;
			break;
		default:
			LOG_FATAL << "1. ProcessType = " << r.which() << " not exists";
			exit(-1);
			break;
		}
		P_total += dR * dt_cell;
		channel++;
	}
	// P_total is the total rate of reation within dt
	if (P_total > 0.15 && !type2_warned)
	{
		LOG_WARNING << "P(dt) = " << P_total << " may be too large";
		type2_warned = true;
	}
	if (Srandom::init_dis(Srandom::gen) > P_total)
	{
		channel = -1;
	}
	else
	{
		// Normalize P_channels using P_total to sample
		for (auto &item : P_channels)
			item /= P_total;
		double p = Srandom::init_dis(Srandom::gen);
		for (int i = 0; i < P_channels.size(); ++i)
		{
			if (P_channels[i] > p)
			{
				channel = i;
				break;
			}
		}
	}

	//LOG_INFO<<"Lido channel: "<<channel;

	// If a scattering happens:
	if (channel >= 0)
	{
		if (channel >= AllProcesses[absid].size())
		{
			LOG_INFO << p00 <<"//"<< pIn.p << "//" << temp;
			LOG_FATAL << "3. Channel = " << channel << " not exists";
			exit(-1);
		}
		// Final state holder FS
		std::vector<fourvec> FS;
		// Sampe its differential final state
		switch (AllProcesses[absid][channel].which())
		{
		case 0:
			boost::get<Rate22>(AllProcesses[absid][channel]).sample({E_cell, temp}, FS);
			break;
		case 1:
			boost::get<Rate23>(AllProcesses[absid][channel]).sample({E_cell, temp}, FS);
			break;
		case 2:
			boost::get<Rate32>(AllProcesses[absid][channel]).sample({E_cell, temp}, FS);
			break;
		case 3:
			boost::get<Rate12>(AllProcesses[absid][channel]).sample({E_cell, temp}, FS);
			break;
		case 4:
			boost::get<Rate21>(AllProcesses[absid][channel]).sample({E_cell, temp}, FS);
			break;
		default:
			LOG_FATAL << "2. Channel = " << channel << " not exists";
			exit(-1);
			break;
		}
		// rotate the final state back and boost it back to the lab frame
		for (auto &pmu : FS)
		{
			pmu = pmu.rotate_back(p_cell);
			pmu = pmu.boost_back(v3cell[0], v3cell[1], v3cell[2]);
		}

		// elastic process changes the momentum immediately
		if (channel == 0 || channel == 1)
		{
			pIn.p = FS[0];
			particle ep;
			if (FS[1].boost_to(v3cell[0], v3cell[1], v3cell[2]).t() >  Lido_Ecut){
				double tempid;
				if (channel == 1) tempid = 21;
				else tempid = Srandom::sample_flavor(3);
				particle ep = produce_parton(tempid, pIn, FS[1], pIn.x, temp, v3cell, false, true);
				pOut_list.push_back(ep);
			}
				
		}
		// inelastic process takes a finite time to happen
		if (channel == 2 || channel == 3)
		{
			// add the radiating daughters as ``virtual" particles
			// these are virtual particles, that pops out at a rate
			// given by the incoherent calculation, which will be
			// dropped (suppressed) according to the LPM effect
			particle vp = produce_parton(21, pIn, FS[2], x0, temp, v3cell);

			// The local 2->2 mean-free-path is estimated with
			// the qhat_hard integrate from the 2->2 rate
			double xfrac = vp.p.t() / pIn.p.t();
			double local_qhat = 0.;
			double vp_cell = vp.p.boost_to(v3cell[0], v3cell[1], v3cell[2]).t();
			double boost_factor = vp.p.t() / vp_cell;
			BOOST_FOREACH (Process &r, AllProcesses[21])
			{
				switch (r.which())
				{
				case 0:
					if (boost::get<Rate22>(r).IsActive())
					{
						tensor A = boost::get<Rate22>(r).GetSecondM(
							{vp_cell, temp});
						local_qhat += A.T[1][1] + A.T[2][2];
					}
					break;
				default:
					break;
				}
			}
			double mD2 = t_channel_mD2->get_mD2(temp);
			// estimate mfp in the lab frame
			vp.mfp0 = LPM_prefactor * mD2 / local_qhat * boost_factor;
			pIn.radlist.push_back(vp);
			if (FS[1].boost_to(v3cell[0], v3cell[1], v3cell[2]).t() > Lido_Ecut){	
				double tempid;
				if (channel == 3) tempid = 21;
				else tempid = Srandom::sample_flavor(3);
				particle ep = produce_parton(tempid, pIn, FS[1], pIn.x, temp, v3cell, false, true);
				pOut_list.push_back(ep);
			}
				
		}
		if (channel == 4 || channel == 5)
		{
			// Absorption processes happens mostly for gluon energy ~ 3*T, therefore we negelected the LPM effect
			pIn.p = FS[0];
		}
		if (channel == 6)
		{
			// add the radiating daughters as ``virtual" particles
			// these are virtual particles, that pops out at a rate
			// given by the incoherent calculation, which will be
			// dropped (suppressed) according to the LPM effect
			particle vp = produce_parton(21, pIn, FS[1], x0, temp, v3cell);

			double mD2 = t_channel_mD2->get_mD2(temp);
			// estimate mfp in the lab frame
			vp.mfp0 = LPM_prefactor * mD2 / qhatg * pIn.p.t() / E_cell;
			pIn.radlist.push_back(vp);
		}
		if (channel == 7)
		{
			// Absorption processes happens mostly for gluon energy ~ 3*T, therefore we negelected the LPM effect
			pIn.p = FS[0];
		}
		// inelastic process takes a finite time to happen
		if (channel == 8 || channel == 9)
		{
			// This should only happen for gluon
			// gluon splits to q+qbar, pid changing!!
			particle vp = produce_parton(Srandom::sample_flavor(3), pIn, FS[2], x0, temp, v3cell);

			// The local 2->2 mean-free-path is estimated with
			// the qhat_hard integrate from the 2->2 rate
			double xfrac = vp.p.t() / pIn.p.t();
			double local_qhat = 0.;
			double vp_cell = vp.p.boost_to(v3cell[0], v3cell[1], v3cell[2]).t();
			double boost_factor = vp.p.t() / vp_cell;
			BOOST_FOREACH (Process &r, AllProcesses[21])
			{
				switch (r.which())
				{
				case 0:
					if (boost::get<Rate22>(r).IsActive())
					{
						tensor A = boost::get<Rate22>(r).GetSecondM(
							{vp_cell, temp});
						local_qhat += A.T[1][1] + A.T[2][2];
					}
					break;
				default:
					break;
				}
			}
			double mD2 = t_channel_mD2->get_mD2(temp);
			// estimate mfp in the lab frame
			vp.mfp0 = LPM_prefactor * mD2 / local_qhat * boost_factor;
			pIn.radlist.push_back(vp);
		}
		if (channel == 10)
		{
			// This should only happen for gluon
			// gluon splits to q+qbar, pid changing!!
			particle vp = produce_parton(Srandom::sample_flavor(3), pIn, FS[1], x0, temp, v3cell);

			double mD2 = t_channel_mD2->get_mD2(temp);
			// estimate mfp in the lab frame
			vp.mfp0 = LPM_prefactor * mD2 / qhatg * pIn.p.t() / E_cell;
			pIn.radlist.push_back(vp);
		}
	}

	// Handle the virtual particle, (only for real mother parton)
	if ((!pIn.radlist.empty()) && (!pIn.is_virtual))
	{
		// loop over each virtual particles
		for (std::vector<particle>::iterator it = pIn.radlist.begin(); it != pIn.radlist.end();)
		{
			int split_type = 0;
			if (pIn.pid != 21 && it->pid == 21)
				split_type = 1; // q --> q + g
			if (pIn.pid == 21 && it->pid == 21)
				split_type = 2; // g --> g + g
			if (pIn.pid == 21 && it->pid != 21)
				split_type = 3; // g --> q + qbar
				
				
			double taun = formation_time(it->mother_p, it->p, temp, split_type);
			std::vector<particle> pnew_Out;
			double tau = std::sqrt(it->x.t()*it->x.t()-it->x.z()*it->x.z());
			double dt_daughter = calcualte_dt_from_dtau(it->x, it->p, 
					tau, dtau); 
			
			// In the local (adiabatic LPM) mode,
			// use local information to determine the number of rescatterings
			bool Adiabatic_LPM = true;
			if (Adiabatic_LPM){
			    do{
			        update_particle_momentum_Lido(dt_daughter, temp, 
											v3cell, *it, pnew_Out);
					taun = formation_time(it->mother_p, it->p, temp, split_type);
			    }while(it->x.t() - it->x0.t() <= taun);
			    // once finished, reset its x to where it is produced
			    it->x = it->x0;
			}
			
			// update the partilce for this time step
			update_particle_momentum_Lido(dt_daughter, temp, 
											v3cell, *it, pnew_Out);
			
			if (it->x.t() - it->x0.t() <= taun && !Adiabatic_LPM){
				it++;
			}
			else{
				// if formed, apply LPM suppression
				double Acceptance = 0.;
				if (!it->is_vac){
					// for medium-induced radiation
					// 1): a change of running-coupling from elastic broadening
					double kt20 = measure_perp(it->mother_p, it->p0).pabs2();
					double kt2n = measure_perp(it->mother_p, it->p).pabs2();
					double Running = alpha_s(kt2n, it->T0) / alpha_s(kt20, it->T0);
					// 2): a dead-cone approximation for massive particles
					double theta2 = kt2n / std::pow(it->p.t(), 2);
					double thetaM2 = std::pow(pIn.mass / it->mother_p.t(), 2);
					double DeadCone = std::pow(theta2 / (theta2 + thetaM2), 2);
					// 3): an NLL-inspired suppression factor for the LPM effect
					double mD2 = t_channel_mD2->get_mD2(temp);
					double lnQ2_1 = std::log(1. + taun / it->mfp0);
					double lnQ2_0 = std::log(1. + 6 * it->p.t() * temp / mD2);
					double log_factor = std::sqrt(lnQ2_1 / lnQ2_0);
					double LPM = std::min(it->mfp0 / taun * log_factor, 1.);
					// The final rejection factor
					Acceptance = LPM * Running * DeadCone;
				}
				else{
					// for vacuum-like radiaiton formed inside the medium
					double kt20 = measure_perp(it->mother_p, it->p0).pabs2();
					double Deltakt2 = measure_perp(it->mother_p, it->p - it->p0).pabs2();
					if (Deltakt2 > Rvac * kt20)
						Acceptance = 0.0;
					else
						Acceptance = 1.0;
				}
				if (Srandom::rejection(Srandom::gen) < Acceptance){
					// accepted branching causes physical effects
					// momentum change, and put back on shell
					//double xx = it->p0.t()/it->mother_p.t();
					//pIn.p = pIn.p * (1. - xx);
					pIn.p = pIn.p - it->p;
					pIn.p.a[0] = std::sqrt(pIn.p.pabs2() + pIn.mass * pIn.mass);
					// for g -> q + qbar, pid change
					// also discard all other pre-splitting: causes higher-order difference
					if (split_type == 3)
					{
						pIn.pid = -it->pid;
						pIn.radlist.clear();
					}
					// label it as real and put it in output particle list
					it->is_virtual = false;
					ofs<<it->p<<" "<<it->mother_p<<" "<<it->p0<<" "<<it->x<<" "<<it->x0<<" "<<it->pid<<" "<<it->is_vac<<" "<<it->is_recoil<<" "<<it->T0 <<" "<<it->weight<<std::endl;
					pOut_list.push_back(*it);
				}
				// remove it from the radlist
				it = pIn.radlist.erase(it);
			}
		}
	}

	// For virtual particle, we rescale the energy back to the initial value
	// so that we only focus on the elastic broadening effect.
	//if (pIn.is_virtual) pIn.p = pIn.p*(pIn.p0.t()/pIn.p.t());

	// Add the mother parton to the end of the output list
	//also transform back its pid
	pIn.pid = temp_pid;
	pOut_list.push_back(pIn);
	
	return pOut_list.size();
}

particle produce_parton(int pid, particle &mother_parton, fourvec vp0, fourvec vx0, double T, std::vector<double> &v3cell, bool is_virtual, bool is_recoil)
{
	particle vp;
	vp.pid = pid;
	vp.mass = 0.0;
	vp.weight = mother_parton.weight;
	vp.p0 = vp0;
	vp.p = vp.p0;
	vp.x0 = vx0;
	vp.x = vp.x0;
	vp.mother_p = mother_parton.p;
	vp.T0 = T;
	vp.is_vac = false;
	vp.is_virtual = is_virtual;
	vp.is_recoil = is_recoil;
	vp.vcell.resize(3);
	vp.vcell[0] = v3cell[0];
	vp.vcell[1] = v3cell[1];
	vp.vcell[2] = v3cell[2];

	return vp;
}

void output(const std::vector<particle> plist, std::string fname)
{
	int i = 0;
	std::ofstream f(fname);
	for (auto &p : plist)
	{
		f << std::setw(10) << std::setfill(' ') << i << "  "	  // particle index, i10,2x
		  << std::setw(10) << std::setfill(' ') << p.pid << "  "; // particle id, i10,2x
		f << p.p.x() << "  "
		  << p.p.y() << "  "
		  << p.p.z() << "  "
		  << p.p.t() << "  "
		  << p.mass << "  "
		  << p.x.x() / fmc_to_GeV_m1 << "  "
		  << p.x.y() / fmc_to_GeV_m1 << "  "
		  << p.x.z() / fmc_to_GeV_m1 << "  "
		  << p.x.t() / fmc_to_GeV_m1 << " "
		  << p.x0.x() / fmc_to_GeV_m1 << "  "
		  << p.x0.y() / fmc_to_GeV_m1 << "  "
		  << p.x0.z() / fmc_to_GeV_m1 << "  "
		  << p.x0.t() / fmc_to_GeV_m1 << "\n";
		i++;
	}
}

void output_oscar(const std::vector<particle> plist, int abspid, std::string fname)
{
	// output OSCAR Format
	int Nparticles = 0;
	for (auto &p : plist)
	{
		if (std::abs(p.pid) == abspid)
			Nparticles++;
	}

	std::ofstream f(fname);
	f << "OSC1997A\n";
	f << "final_id_p_x\n";
	f << "     lbt  1.0alpha   208     82   208     82   aacm  0.1380E+04         1\n";
	f << "         1  " << std::setw(10)
	  << Nparticles << "     0.001     0.001     1     1        1\n";

	int i = 0;
	for (auto &p : plist)
	{
		if (std::abs(p.pid) == abspid)
		{
			f << std::setw(10) << std::setfill(' ') << i << "  "	  // particle index, i10,2x
			  << std::setw(10) << std::setfill(' ') << p.pid << "  "; // particle id, i10,2x
			f << ff(p.p.x()) << "  "
			  << ff(p.p.y()) << "  "
			  << ff(p.p.z()) << "  "
			  << ff(p.p.t()) << "  "
			  << ff(p.mass) << "  "
			  << ff(p.x.x() / fmc_to_GeV_m1) << "  "
			  << ff(p.x.y() / fmc_to_GeV_m1) << "  "
			  << ff(p.x.z() / fmc_to_GeV_m1) << "  "
			  << ff(p.x.t() / fmc_to_GeV_m1) << "  "
			  << ff(p.Tf) << "  "
			  << ff(p.vcell[0]) << "  "
			  << ff(p.vcell[1]) << "  "
			  << ff(p.vcell[2]) << "  "
			  << ff(p.p0.x()) << "  "
			  << ff(p.p0.y()) << "  "
			  << ff(p.p0.z()) << "  "
			  << ff(p.p0.t()) << "  "
			  << ff(p.weight) << "  "
			  << ff(0.0) << "\n";
		}
		i++;
	}
}

double calcualte_dt_from_dtau(fourvec x, fourvec p, double tau, double dtau)
{
	double vz = p.z() / p.t();
	double t_m_zvz = x.t() - x.z() * vz;
	double one_m_vz2 = 1. - vz * vz;
	double dtau2 = dtau * (dtau + 2. * tau);
	double dt_lab = (std::sqrt(t_m_zvz * t_m_zvz + one_m_vz2 * dtau2) - t_m_zvz) / one_m_vz2;
	return dt_lab;
}

double mean_pT(const std::vector<particle> plist)
{
	double W = 0., wsum_pT = 0.;
	for (auto &p : plist)
	{
		wsum_pT += p.weight * std::sqrt(p.p.x() * p.p.x() + p.p.y() * p.p.y());
		W += p.weight;
	}
	return wsum_pT / W;
}
double mean_E(const std::vector<particle> plist)
{
	double W = 0., wsum_pT = 0.;
	for (auto &p : plist)
	{
		wsum_pT += p.weight * p.p.t();
		W += p.weight;
	}
	return wsum_pT / W;
}
