#include <string>
#include <iostream>
#include <exception>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <fstream>

#include "TableBase.h"
#include "Xsection.h"
#include "Rate.h"
#include "matrix_elements.h"

void test_r(void);
void test_x(void);
void test_config(void);
void test_table(void);

int main(){
	test_r();
	//test_x();
	//test_config();
	//test_table();
	return 0;
}
void test_r(void){
	using boost::property_tree::ptree;
    ptree config;
    std::ifstream input("settings.xml");
    read_xml(input, config);
	double mu = config.get_child("Boltzmann.QCD").get<double>("mu");
	std:: cout <<  "mu = " << mu << std::endl;
    initialize_mD_and_scale(0, mu);

	auto rQq2Qq = std::make_shared<Rate<2, 2, double(*)(double, void*)>>
	("Qq2Qq", config.get_child("Boltzmann"), dX_Qq2Qq_dt);
	rQq2Qq->init();
	rQq2Qq->save("qhat.h5");

	auto rQg2Qg = std::make_shared<Rate<2, 2, double(*)(double, void*)>>
	("Qg2Qg", config.get_child("Boltzmann"), dX_Qg2Qg_dt);
	rQg2Qg->init();
	rQg2Qg->save("qhat.h5");

	//auto r = std::make_shared<Rate<2, 2, double(*)(double, void*)>>
	//("resonance", config.get_child("Boltzmann"), dX_res_dt);
	//r->init();
	//r->save("table.h5");
}
/*
void test_x(void){
	using boost::property_tree::ptree;
    ptree config;
    std::ifstream input("settings.xml");
    read_xml(input, config);
	auto xQq2Qq = std::make_shared<Xsection<2, double(*)(double, void*)>>
	("Qq2Qq/xsection", config.get_child("Boltzmann"), dX_Qq2Qq_dt);
	auto xQg2Qg = std::make_shared<Xsection<2, double(*)(double, void*)>>
	("Qg2Qg/xsection", config.get_child("Boltzmann"), dX_Qg2Qg_dt);
	xQq2Qq->init();
	xQq2Qq->save("table.h5");
	xQg2Qg->init();
	xQg2Qg->save("table.h5");
}
void test_config(void){
    using boost::property_tree::ptree;
    ptree pt;
    std::ifstream input("settings.xml");
    read_xml(input, pt);
    for(ptree::value_type const& v: pt.get_child("Boltzmann") ) {
		std::cout << std::endl << v.first << std::endl;
        if( v.first > "process" ) {
           std::cout << v.second.get<std::string>("name") << " "
					<< v.second.get("<xmlattr>.status", "-") 
					<< std::endl
         			<< v.second.get<std::string>("probe") << std::endl;
			std::cout << "x-table info:" << std::endl;
			for(auto & vv: pt.get_child("Boltzmann."+v.first+".xsection") ){
				std::cout << "\t" << vv.first << " " 
						   << vv.second.get("<xmlattr>.name", "none") << " "
							<< vv.second.data() << " "
							<< std::endl;
			}
			std::cout << "r-table info:" << std::endl;;
			for(auto & vv: pt.get_child("Boltzmann."+v.first+".rate")){
				std::cout << "\t" << vv.first << " " 
						   << vv.second.get("<xmlattr>.name", "none") << " "
							<< vv.second.data() << " "
							<< std::endl;
			}
		}
    }
}

void test_table(void){
	TableBase<scalar,2> T1(std::string("Table1"), {{5,5}}, {{0,0}}, {{1,1}});
    for (auto i=0; i<5; ++i) {
   		for (auto j=0; j<5; ++j) {
   			T1.SetTableValue({i, j}, {i*j});
    	}
    }
    std::cout << T1.InterpolateTable({.25,.35}) << std::endl;

    TableBase<scalar,3> T2(std::string("Table2"), {{5,5,10}}, {{0,0,0}}, {{1,1,3}});
    for (auto i=0; i<5; ++i) {
   		for (auto j=0; j<5; ++j) {
   		    for (auto k=0; k<10; ++k) {
	   			T2.SetTableValue({i, j, k}, {i*j/(k+1.)});
			}
    	}
    }
    std::cout << T2.InterpolateTable({.3,.3,.0}) << std::endl;
    
    TableBase<fourvec,3> T3(std::string("Table3"), {{5,5,10}}, {{0,0,0}}, {{1,1,3}});
    for (auto i=0; i<5; ++i) {
   		for (auto j=0; j<5; ++j) {
   		    for (auto k=0; k<10; ++k) {
	   			T3.SetTableValue({i, j, k}, {0,i,j,k});
			}
    	}
    }
    std::cout << T3.InterpolateTable({.3,.3,.3}).boost_to(0.99, 0., 0.) << std::endl;
    
    
    T3.Save("field.h5");
    TableBase<fourvec,3> T4(std::string("Table3"), {{5,5,10}}, {{0,0,0}}, {{1,1,3}});    
    T4.Load("field.h5");
    T4.Save("field-2.h5");

}*/

