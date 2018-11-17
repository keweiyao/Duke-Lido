#include <string>
#include <iostream>
#include <exception>

#include "simpleLogger.h"
#include "workflow.h"

// This sample program evolve heavy quark (E0=30, M=1.3) in a static medium for t=3.0fm/c and calculate the average energy loss (<E-E0>)
int main(int argc, char* argv[]){
	if (argc < 7){
                std::cout << "Usage: Lido-TabGen <lido-setting> <tab-path> <mu> <afix> <A> <B>" << std::endl;
                exit(-1);
        }
        std::string lido_setting(argv[1]);
	std::string table_path(argv[2]);
        double mu = atof(argv[3]);
        double afix = atof(argv[4]);
        double A = atof(argv[5]);
        double B = atof(argv[6]);

	initialize("new", lido_setting, table_path, mu, afix, A, B);
	return 0;
}


