#include "TableBase.h"
#include <string>
#include <iostream>


double f(std::vector<double> X);
fourvec f1(std::vector<double> X);


int main(){
    TableBase<double,2> T1(std::string("Table1"), &f, {{5,5}}, {{0,0}}, {{1,1}});
    for (auto i=0; i<5; ++i) {
   		for (auto j=0; j<5; ++j) {
   			T1.SetTableValue({i, j}, i*j);
    	}
    }
    std::cout << T1.InterpolateTable({.25,.35}) << std::endl;

    TableBase<double,3> T2(std::string("Table2"), &f, {{5,5,10}}, {{0,0,0}}, {{1,1,3}});
    for (auto i=0; i<5; ++i) {
   		for (auto j=0; j<5; ++j) {
   		    for (auto k=0; k<10; ++k) {
	   			T2.SetTableValue({i, j, k}, i*j/(k+1.));
			}
    	}
    }
    std::cout << T2.InterpolateTable({.3,.3,.0}) << std::endl;
    
    TableBase<fourvec,3> T3(std::string("Table3"), &f1, {{5,5,10}}, {{0,0,0}}, {{1,1,3}});
    for (auto i=0; i<5; ++i) {
   		for (auto j=0; j<5; ++j) {
   		    for (auto k=0; k<10; ++k) {
	   			T3.SetTableValue({i, j, k}, {0,i,j,k});
			}
    	}
    }
    std::cout << T3.InterpolateTable({.3,.3,.3}).boost_to(0.99, 0., 0.) << std::endl;
}

double f(std::vector<double> X){
  return 1.0;
}

fourvec f1(std::vector<double> X){
	return fourvec{0,1,2,3.5};
}
