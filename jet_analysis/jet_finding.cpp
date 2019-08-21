#include <fstream>
#include <string>
#include <math.h>
#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/FastJet3.h"
#include "fastjet/Selector.hh"
using namespace Pythia8;

template <typename T>
std::string to_string(T value)
{
	std::ostringstream os;
	os << value;
	return os.str();
}

template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 2)
{
	std::ostringstream out;
	out << std::setprecision(n) << a_value;
	return out.str();
}


void simu(double eCM, double jetR, double ymin, double ymax, int nEvent)
{
	std::fstream fs;
	fs.open("./jet_results/jet_" + to_string(eCM) + "_" + to_string(jetR) + "_" + to_string(ymin) + "_" + to_string(ymax) + ".dat", std::fstream::out);
	std::fstream fs2;
	fs2.open("./jet_results/bjet_" + to_string(eCM) + "_" + to_string(jetR) + "_" + to_string(ymin) + "_" + to_string(ymax) + ".dat", std::fstream::out);
	std::fstream fs3;
	fs3.open("./jet_results/cjet_" + to_string(eCM) + "_" + to_string(jetR) + "_" + to_string(ymin) + "_" + to_string(ymax) + ".dat", std::fstream::out);

    int power = -1; // power=-1: anti-kT
	double jetRadius = jetR, jetpTMin = 10.;
	double jetyMin = ymin, jetyMax = ymax;

	std::vector<double> pTBin({30, 40, 50, 56, 63, 70, 79, 89, 100, 112, 125, 141, 158, 177, 199, 251, 398});
	//std::vector<double> pTBin({50,100,200,300,400});
	std::vector<double> jet_cs(pTBin.size(), 0.);
	std::vector<double> err(pTBin.size(), 0.);
	std::vector<double> sqSum(pTBin.size(), 0.);
	std::vector<double> bjet_cs(pTBin.size(), 0.);
	std::vector<double> berr(pTBin.size(), 0.);
	std::vector<double> bsqSum(pTBin.size(), 0.);
	std::vector<double> cjet_cs(pTBin.size(), 0.);
	std::vector<double> cerr(pTBin.size(), 0.);
	std::vector<double> csqSum(pTBin.size(), 0.);

	vector<double> pTHatBin({10, 20, 30, 50, 70, 90, 110, 150, 200, 300, 500, 700});
	size_t nBin = pTHatBin.size() - 1;

	//fastjet.init();
	fastjet::JetDefinition jetDef(fastjet::genkt_algorithm, jetRadius, power);
	std::vector<fastjet::PseudoJet> fjInputs;
	fastjet::Selector select_rapidity = fastjet::SelectorRapRange(jetyMin, jetyMax);

	for (size_t iBin = 0; iBin < nBin; iBin++)
	{
		std::vector<int> jet_ct = std::vector<int>(pTBin.size() - 1, 0);
		std::vector<int> bjet_ct = std::vector<int>(pTBin.size() - 1, 0);
		std::vector<int> cjet_ct = std::vector<int>(pTBin.size() - 1, 0);

		double sigmaGen=0.;

		for (int iEvent = 0; iEvent < nEvent; ++iEvent)
		{

			fjInputs.resize(0);

			std::string file_name = "./results/" + to_string_with_precision(pTHatBin[iBin], 3) + "-" + to_string_with_precision(pTHatBin[iBin + 1], 3) + "-" + to_string_with_precision(iEvent, 3);

			std::cout<< file_name<<std::endl;

			std::ifstream infile(file_name);
			if (infile.is_open())
			{
				int pid;
				double px;
				double py;
				double pz;
				double e;
				std::string str;
				char c;
				while (infile >> pid >> e >> px >> py >> pz >> sigmaGen)
				{
					fjInputs.push_back(fastjet::PseudoJet(px, py, pz, e));
					//determine channel for final heavy quarks
					fjInputs[fjInputs.size() - 1].set_user_index(pid);
				}
			}
			infile.close();

			/* 
			double total_pT=0., hard_pT=0.;

			std::vector<fastjet::PseudoJet> temp;
			for (int iter=0; iter<fjInputs.size(); iter++)
			{
				total_pT+=fjInputs[iter].pt();
				if (fjInputs[iter].pt()>1)
				{
					temp.push_back(fjInputs[iter]);
					hard_pT+=fjInputs[iter].pt();
				}
			}

			fjInputs=temp;

			cout<<"total pT: "<<total_pT <<", hard pT: "<<hard_pT<<endl;
			*/

			//cout<<"parton number: "<<fjInputs.size()<<endl;
			vector<fastjet::PseudoJet> inclusiveJets, jets;
			fastjet::ClusterSequence clustSeq(fjInputs, jetDef);
			inclusiveJets = clustSeq.inclusive_jets(jetpTMin);
			std::vector<fastjet::PseudoJet> selected_jets = select_rapidity(inclusiveJets);
			jets = sorted_by_pt(selected_jets);
			for (int k = 0; k < jets.size(); k++)
			{
				vector<fastjet::PseudoJet> constituents = jets[k].constituents();

				for (int p = 0; p < constituents.size(); p++)
				{
					std::string pid = to_string(constituents[p].user_index());
					if (pid.find("5") != std::string::npos)
					{
						//cout<<pid<<endl;
						for (unsigned int j = 0; j < pTBin.size() - 1; j++)
						{
							if (jets[k].pt() >= pTBin[j] && jets[k].pt() < pTBin[j + 1])
							{
								bjet_ct[j]++;
								break;
							}
						}
						break;
					}
					if (pid.find("4") != std::string::npos)
					{
						//cout<<pid<<endl;
						for (unsigned int j = 0; j < pTBin.size() - 1; j++)
						{
							if (jets[k].pt() >= pTBin[j] && jets[k].pt() < pTBin[j + 1])
							{
								cjet_ct[j]++;
								break;
							}
						}
						break;
					}
				}
				for (unsigned int j = 0; j < pTBin.size() - 1; j++)
				{
					if (jets[k].pt() >= pTBin[j] && jets[k].pt() < pTBin[j + 1])
					{
						jet_ct[j]++;
						break;
					}
				}
			}
		}

		double sigmapb_weight = sigmaGen * 1.0e9 / nEvent;

		for (unsigned int j = 0; j < pTBin.size() - 1; j++)
		{
			jet_cs[j] += jet_ct[j] * sigmapb_weight;
			sqSum[j] += jet_ct[j] * pow(sigmapb_weight, 2);
			bjet_cs[j] += bjet_ct[j] * sigmapb_weight;
			bsqSum[j] += bjet_ct[j] * pow(sigmapb_weight, 2);
			cjet_cs[j] += cjet_ct[j] * sigmapb_weight;
			csqSum[j] += cjet_ct[j] * pow(sigmapb_weight, 2);
		}
	}
	for (unsigned int j = 0; j < pTBin.size() - 1; j++)
	{
		err[j] = jet_cs[j] / sqrt(pow(jet_cs[j], 2) / sqSum[j]);
		cerr[j] = cjet_cs[j] / sqrt(pow(cjet_cs[j], 2) / csqSum[j]);
		berr[j] = bjet_cs[j] / sqrt(pow(bjet_cs[j], 2) / bsqSum[j]);
		cout << (pTBin[j] + pTBin[j + 1]) / 2 << " " << jet_cs[j] / (pTBin[j + 1] - pTBin[j]) << " " << err[j] / (pTBin[j + 1] - pTBin[j]) << endl;
		fs << (pTBin[j] + pTBin[j + 1]) / 2 << " " << jet_cs[j] / (pTBin[j + 1] - pTBin[j]) << " " << err[j] / (pTBin[j + 1] - pTBin[j]) << endl;
		fs2 << (pTBin[j] + pTBin[j + 1]) / 2 << " " << bjet_cs[j] / (pTBin[j + 1] - pTBin[j]) << " " << berr[j] / (pTBin[j + 1] - pTBin[j]) << endl;
		fs3 << (pTBin[j] + pTBin[j + 1]) / 2 << " " << cjet_cs[j] / (pTBin[j + 1] - pTBin[j]) << " " << cerr[j] / (pTBin[j + 1] - pTBin[j]) << endl;
	}

	fs.close();
	fs2.close();
	fs3.close();
}

int main(int argc, char **argv)
{
	simu(std::stod(argv[1]), std::stod(argv[2]), std::stod(argv[3]), std::stod(argv[4]),std::stoi(argv[5]));
	//simu(5020,0.4,0,2.8);

	return 0;
}
