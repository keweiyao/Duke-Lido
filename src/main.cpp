#include "TableBase.h"
#include <string>
#include <iostream>


double f(std::vector<double> X);
VectorT f1(std::vector<double> X);


int main(){
    TableBase<double,2> T1(std::string("Table1"), &f, {{2,2}});
    T1.SetTableValue({0, 0}, 1.0);
	T1.SetTableValue({0, 1}, 2.0);
	T1.SetTableValue({1, 0}, 3.0);
	T1.SetTableValue({1, 1}, 4.0);
    std::cout << T1.InterpolateTable({.3,.3}) << std::endl;

    TableBase<double,3> T2(std::string("Table2"), &f, {{2,2,2}});
    T2.SetTableValue({0, 0, 0}, 1.0);
	T2.SetTableValue({0, 1, 0}, 2.0);
	T2.SetTableValue({1, 0, 0}, 3.0);
	T2.SetTableValue({1, 1, 0}, 4.0);
    T2.SetTableValue({0, 0, 1}, 5.0);
	T2.SetTableValue({0, 1, 1}, 6.0);
	T2.SetTableValue({1, 0,1}, 7.0);
	T2.SetTableValue({1, 1, 1}, 8.0);
    std::cout << T2.InterpolateTable({.3,.3,.5}) << std::endl;
}

double f(std::vector<double> X){
  return 1.0;
}

VectorT f1(std::vector<double> X){
	return VectorT{0,1,2,3.5};
}
