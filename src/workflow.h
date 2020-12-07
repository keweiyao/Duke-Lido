#ifndef WORKFLOW_H
#define WORKFLOW_H

#include <vector>



void initialize(std::string mode, std::string setting_path, std::string table_path,
	double mu, double afix, double theta, double cut);

int update_particle_momentum_Lido(double dt, double temp, std::vector<double> v3cell, particle & pIn, std::vector<particle> & pOut_list, double Tf);


void output_oscar(const std::vector<particle> plist, int abspid, std::string fname);
void SampleFlavorAndColor(int mother_id, int mother_col, int mother_acol, 
                int channel, int daughter_id, 
                int & daughter_col, int & daughter_acol,
                int & new_mother_col, int & new_mother_acol);
#endif
