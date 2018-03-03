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
#include "simpleLogger.h"

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
	std::vector<fourvec> FS;
	using boost::property_tree::ptree;
    ptree config;
    std::ifstream input("settings.xml");
    read_xml(input, config);
	double mu = config.get_child("Boltzmann.QCD").get<double>("mu");
    initialize_mD_and_scale(0, mu);

	auto rQq2Qq = std::make_shared<Rate<2, 2, double(*)(const double, void*)>>
	("Boltzmann/Qq2Qq", "settings.xml", dX_Qq2Qq_dt);
	rQq2Qq->init();
	rQq2Qq->save("table.h5");

	auto rQq2Qqg = std::make_shared<Rate<3, 3, double(*)(const double*, void*)>>
	("Boltzmann/Qq2Qqg", "settings.xml", M2_Qq2Qqg);
	rQq2Qqg->init();
	LOG_INFO << "-----Done-----";
	rQq2Qqg->save("table.h5");
    
	for(int i=0; i<100000; i++){
		if (i%10000==0) LOG_INFO << i << " scattering events sampled";
		double T = 0.15 + std::rand()*0.85/RAND_MAX;
		double E = 5. + std::rand()*25./RAND_MAX;
		double dt = 1.0 + std::rand()*5.0/RAND_MAX;
		rQq2Qqg->sample({E,T,dt},FS);
		rQq2Qq->sample({E,T},FS);
	}

	/*
	std::vector<fourvec> plist(10000);
	double T = 0.6;
	double dt = 0.2;
	for (auto &p : plist) p = fourvec{3, 0., 0., std::sqrt(3*3-1.3*1.3)};
	for (int i=0; i<400; i++){
		for (auto &p : plist){
			double E = p.t();
			double prob = dt*rQq2Qq->GetZeroM({E, T}).s;
			if (std::rand()*1./RAND_MAX < prob){
				rQq2Qq->sample({E, T}, FS);
				p = FS[0].rotate_back(p);
			}
		}
		if(i%20 ==0)
			for (auto &p : plist) std::cout << p << std::endl;
	}*/
 

	//auto rQg2Qg = std::make_shared<Rate<2, 2, double(*)(double, void*)>>
	//("Qg2Qg", config.get_child("Boltzmann"), dX_Qg2Qg_dt);
	//rQg2Qg->init();
	//rQg2Qg->save("table.h5");

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
}*/

scalar approx(std::vector<double> x){
	return scalar{x[0]*x[1]*x[1]};
}
void test_table(void){
	TableBase<scalar,2> T1(std::string("Table1"), {{4,4}}, {{0,0}}, {{1,1}});
    for (auto i=0; i<4; ++i) {
   		for (auto j=0; j<4; ++j) {
   			T1.SetTableValue({i, j}, {i*j*j});
   			std::cout << i*j*j << " ";
    	}
    	std::cout << std::endl;
    }
    T1.SetApproximateFunction(approx);
    std::cout << T1.InterpolateTable({.5,.5}) << std::endl;

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

}

