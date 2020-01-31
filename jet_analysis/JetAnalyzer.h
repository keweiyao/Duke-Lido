#include <fstream>
#include <string>
#include <math.h>
#include <sstream>
#include <algorithm>
#include <iterator>
#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/FastJet3.h"
#include "fastjet/Selector.hh"
using namespace Pythia8;

void splitToDouble(const string str, vector<double> &rst, char delim = ' ')
{
	std::stringstream ss(str);
	std::istream_iterator<std::string> begin(ss);
	std::istream_iterator<std::string> end;
	std::vector<std::string> vstrings(begin, end);

	for (auto &str : vstrings)
	{
		rst.push_back(stod(str));
	}
}

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

class JetAnalyzer{

    protected:
    std::string analyzeType;
    std::fstream fs;
    std::vector<double> outputBins;
    std::vector<double> rst;
    std::vector<double> sqSum;
    std::vector<double> err;
    int clusterPower;
    double jetRadius;
    double jetpTMin;
    double jetyMin;
    double jetyMax;

    public: 
    JetAnalyzer(std::string analyzeType, std::vector<double> outputBins, double clusterPower, double jetRadius
   , double jetpTMin, double jetyMin, double jetyMax): analyzeType(analyzeType), outputBins(outputBins), 
   clusterPower(clusterPower), jetRadius(jetRadius), jetpTMin(jetpTMin), jetyMin(jetyMin), jetyMax(jetyMax) {
fs.open("./jet_results/" + analyzeType + "_p=" + to_string(clusterPower)+ "_R=" + to_string(jetRadius) + "_y=" + to_string(jetyMin) + "-" + to_string(jetyMax) + "_pTMin" + to_string(jetpTMin) + ".dat", std::fstream::out);

    rst=std::vector<double>(outputBins.size(), 0.);
	err=std::vector<double>(outputBins.size(), 0.);
    sqSum=std::vector<double>(outputBins.size(), 0.);
   };

   void doJetFinding(std::vector<fastjet::PseudoJet> fjInputs) {}

   void outputResults() {
    for (unsigned int j = 0; j < outputBins.size() - 1; j++)
	{
		err[j] = rst[j] / sqrt(pow(rst[j], 2) / sqSum[j]);
		fs << (outputBins[j] + outputBins[j + 1]) / 2 << " " << rst[j] / (outputBins[j + 1] - outputBins[j]) << " " << err[j] / (outputBins[j + 1] - outputBins[j]) << endl;
	}
    fs.close();
   }

};