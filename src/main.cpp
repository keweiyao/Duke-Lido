#include "TableBase.h"
#include <string>
#include <iostream>

VectorData f(std::vector<double> X);

int main(){
    VectorTable2D T1(std::string("Table1"), &f);
    T1.InterpolateTable(2, 0.3, 0.5);
    T1.SetTableAt(2, 10, 11);
    T1.GetTableAt(2, 3, 1);
}

VectorData f(std::vector<double> X){
  double x = X[0];
  double y = X[1];
  VectorData A{x, x+1, y*3, y*3+1};
  return A;
}
