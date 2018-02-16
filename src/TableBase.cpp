#include "TableBase.h"
#include <iostream>

template <typename T, size_t N>
TableBase<T, N>::TableBase(std::string Name, T (*f)(std::vector<double>)):
_Name(Name), _dims(N), _approximating_function(f)
{
	std::cout<<_Name << " dim=" << _dims <<std::endl;
    std::vector<double> x = {1,1};
    std::cout<<"Test function f(1,1) = " << _approximating_function(x) << std::endl;   
}

VectorTable2D::VectorTable2D(std::string Name, VectorData (*f)(std::vector<double>)):
TableBase<VectorData, 2>(Name, f)
{
}

VectorData VectorTable2D::InterpolateTable(size_t n, ...){
   va_list args;
   va_start(args, n);
   std::cout<<"interp ";
   for (size_t i=0; i<n; ++i) {
       auto value = va_arg(args, double);
       std::cout << value << " ";
   }
   std::cout<<std::endl;
   va_end(args); 
}

void VectorTable2D::SetTableAt(size_t n, ...){
   va_list args;
   va_start(args, n);
   std::cout<<"set ";
   for (size_t i=0; i<n; ++i) {
       auto value = va_arg(args, size_t);
       std::cout << value << " ";
   }
   std::cout<<std::endl;
   va_end(args);
}

VectorData VectorTable2D::GetTableAt(size_t n, ...){
   va_list args;
   va_start(args, n);
   std::cout<<"get ";
   for (size_t i=0; i<n; ++i) {
       auto value = va_arg(args, size_t);
       std::cout << value << " ";
   }
   std::cout<<std::endl;
   va_end(args);
}

