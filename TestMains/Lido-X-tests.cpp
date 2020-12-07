#include <string>
#include <iostream>
#include <exception>

#include "Xsection.h"
#include "simpleLogger.h"
#include "predefine.h"
#include "matrix_elements.h"
#include <fstream>
#include <vector>
#include <iterator>


void test_X22_q(void);
void test_X22_g(void);
void test_X22_hq(void);
void test_X23_q(void);
void test_X23_g(void);
void test_X23_hq(void);

int main(){
    initialize_mD_and_scale(0, 2, 0.3, 4., 4.);
    LOG_INFO << "X22-q";
    test_X22_q();
    LOG_INFO << "X22-g";
    test_X22_g();
    LOG_INFO << "X22-hq";
    test_X22_hq();
    LOG_INFO << "X23-q";
    test_X23_q();
    LOG_INFO << "X23-g";
    test_X23_g();
    LOG_INFO << "X23-hq";
    test_X23_hq();
}

void test_X22_q(void){
    Xsection<HS2PP, 2, double(*)(const double, void*)>
    x1("Boltzmann/qg2qg", "lido_setting.xml", dX_qg2qg),
    x2("Boltzmann/qq2qq", "lido_setting.xml", dX_qq2qq),
    x3("Boltzmann/qqbar2qqbar", "lido_setting.xml", dX_qqbar2qqbar_diff),
    x4("Boltzmann/qqbar2ccbar", "lido_setting.xml", dX_qqbar2qqbar_diff),
    x5("Boltzmann/qqbar2bbbar", "lido_setting.xml", dX_qqbar2qqbar_diff);
    std::vector<Xsection<HS2PP, 2, double(*)(const double, void*)> >
    X{x1, x2, x3, x4, x5};
    std::vector<std::string> S{"qg2qg", "qq2qq", "qqbar2qqbar",
                               "qqbar2ccbar", "qqbar2bbbar"};
    for (auto & x : X) x.load("table.h5");
    double sqrts = 10;
    double T = 0.3;
    double lnsqrts = std::log(sqrts);
    int pid=1;
    for (int j=0; j<5; j++){
        std::ofstream f("Xsamples/"+S[j]+".dat");    
        for (int i=0; i<10000; i++){
            std::vector<fourvec> FS;
            std::vector<int> ID;
            X[j].sample({lnsqrts, T}, pid, FS, ID);
            if (FS.size()==2)
                f << ID[0] << " " << FS[0] << " "
                  << ID[1] << " " << FS[1] << std::endl; 
        }
        f.close();
    }
}

void test_X22_g(void){
    Xsection<HS2PP, 2, double(*)(const double, void*)>
    x1("Boltzmann/gg2gg", "lido_setting.xml", dX_gg2gg),
    x2("Boltzmann/gq2gq", "lido_setting.xml", dX_gq2gq),
    x3("Boltzmann/gg2qqbar", "lido_setting.xml", dX_gg2qqbar),
    x4("Boltzmann/gg2ccbar", "lido_setting.xml", dX_gg2qqbar),
    x5("Boltzmann/gg2bbbar", "lido_setting.xml", dX_gg2qqbar);
    std::vector<Xsection<HS2PP, 2, double(*)(const double, void*)> >
    X{x1, x2, x3, x4, x5};
    std::vector<std::string> S{"gg2gg", "gq2gq", "gg2qqbar",
                               "gg2ccbar", "gg2bbbar"};
    for (auto & x : X) x.load("table.h5");
    double sqrts = 10;
    double T = 0.3;
    double lnsqrts = std::log(sqrts);
    int pid = 21;
    for (int j=0; j<5; j++){
        std::ofstream f("Xsamples/"+S[j]+".dat");    
        for (int i=0; i<10000; i++){
            std::vector<fourvec> FS;
            std::vector<int> ID;
            X[j].sample({lnsqrts, T}, pid, FS, ID);
            if (FS.size()==2)
                f << ID[0] << " " << FS[0] << " "
                  << ID[1] << " " << FS[1] << std::endl; 
        }
        f.close();
    }
}

void test_X22_hq(void){
    Xsection<HS2PP, 2, double(*)(const double, void*)>
    x1("Boltzmann/cg2cg", "lido_setting.xml", dX_qg2qg),
    x2("Boltzmann/cq2cq", "lido_setting.xml", dX_qq2qq),
    x3("Boltzmann/bg2bg", "lido_setting.xml", dX_qg2qg),
    x4("Boltzmann/bq2bq", "lido_setting.xml", dX_qq2qq);
    std::vector<Xsection<HS2PP, 2, double(*)(const double, void*)> >
    X{x1, x2, x3, x4};
    std::vector<std::string> S{"cg2cg", "cq2cq", "bg2bg", "bq2bq"};
    for (auto & x : X) x.load("table.h5");
    double sqrts = 10;
    double T = 0.3;
    double lnsqrts = std::log(sqrts);
    int pid;
    for (int j=0; j<4; j++){
        if (j<2) pid=4;
        else pid=5;
        std::ofstream f("Xsamples/"+S[j]+".dat");    
        for (int i=0; i<10000; i++){
            std::vector<fourvec> FS;
            std::vector<int> ID;
            X[j].sample({lnsqrts, T}, pid, FS, ID);
            if (FS.size()==2)
                f << ID[0] << " " << FS[0] << " "
                  << ID[1] << " " << FS[1] << std::endl; 
        }
        f.close();
    }
}


void test_X23_q(void){
    Xsection<HS2PPP, 2, double(*)(const double*, void*)>
    x1("Boltzmann/qg2qgg", "lido_setting.xml", dX_qg2qgg),
    x2("Boltzmann/qq2qqg", "lido_setting.xml", dX_qq2qqg),
    x3("Boltzmann/qg2qqqbar", "lido_setting.xml", dX_qg2qqqbar),
    x4("Boltzmann/qg2qccbar", "lido_setting.xml", dX_qg2qqqbar),
    x5("Boltzmann/qg2qbbbar", "lido_setting.xml", dX_qg2qqqbar);
    std::vector<Xsection<HS2PPP, 2, double(*)(const double*, void*)> >
    X{x1, x2, x3, x4, x5};
    std::vector<std::string> S{"qg2qgg", "qq2qqg", "qg2qqqbar",
                               "qg2qccbar", "qg2qbbbar"};
    for (auto & x : X) x.load("table.h5");
    double sqrts = 10;
    double T = 0.3;
    double lnsqrts = std::log(sqrts);
    int pid=1;
    for (int j=0; j<5; j++){
        std::ofstream f("Xsamples/"+S[j]+".dat");    
        for (int i=0; i<10000; i++){
            std::vector<fourvec> FS;
            std::vector<int> ID;
            X[j].sample({lnsqrts, T}, pid, FS, ID);
            if (FS.size()==3)
                f << ID[0] << " " << FS[0] << " "
                  << ID[1] << " " << FS[1] << " "
                  << ID[2] << " " << FS[2] << std::endl; 
        }
        f.close();
    }
}

void test_X23_g(void){
    Xsection<HS2PPP, 2, double(*)(const double*, void*)>
    x1("Boltzmann/gg2ggg", "lido_setting.xml", dX_gg2ggg),
    x2("Boltzmann/gg2qgqbar", "lido_setting.xml", dX_gg2qgqbar),
    x3("Boltzmann/gg2cgcbar", "lido_setting.xml", dX_gg2qgqbar),
    x4("Boltzmann/gg2bgbbar", "lido_setting.xml", dX_gg2qgqbar),
    x5("Boltzmann/gq2gqg", "lido_setting.xml", dX_gq2gqg),
    x6("Boltzmann/gq2qqqbar", "lido_setting.xml", dX_gq2qqqbar),
    x7("Boltzmann/gq2cqcbar", "lido_setting.xml", dX_gq2qqqbar),
    x8("Boltzmann/gq2bqbbar", "lido_setting.xml", dX_gq2qqqbar);

    std::vector<Xsection<HS2PPP, 2, double(*)(const double*, void*)> >
    X{x1, x2, x3, x4, x5, x6, x7, x8};
    std::vector<std::string> S{"gg2ggg", "gg2qgqbar", "gg2cgcbar",
                               "gg2bgbbar", "gq2gqg", "gq2qqqbar",
                               "gq2cqcbarr", "gq2bqbbar"};
    for (auto & x : X) x.load("table.h5");
    double sqrts = 10;
    double T = 0.3;
    double lnsqrts = std::log(sqrts);
    int pid=21;
    for (int j=0; j<8; j++){
        std::ofstream f("Xsamples/"+S[j]+".dat");    
        for (int i=0; i<10000; i++){
            std::vector<fourvec> FS;
            std::vector<int> ID;
            X[j].sample({lnsqrts, T}, pid, FS, ID);
            if (FS.size()==3)
                f << ID[0] << " " << FS[0] << " "
                  << ID[1] << " " << FS[1] << " "
                  << ID[2] << " " << FS[2] << std::endl; 
        }
        f.close();
    }
}


void test_X23_hq(void){
    Xsection<HS2PPP, 2, double(*)(const double*, void*)>
    x1("Boltzmann/cg2cgg", "lido_setting.xml", dX_qg2qgg),
    x2("Boltzmann/cq2cqg", "lido_setting.xml", dX_qq2qqg),
    x3("Boltzmann/bg2bgg", "lido_setting.xml", dX_qg2qgg),
    x4("Boltzmann/bq2bqg", "lido_setting.xml", dX_qq2qqg);

    std::vector<Xsection<HS2PPP, 2, double(*)(const double*, void*)> >
    X{x1, x2, x3, x4};
    std::vector<std::string> S{"cg2cgg", "cq2cqg", "bg2bgg", "bq2bqg"};
    for (auto & x : X) x.load("table.h5");
    double sqrts = 10;
    double T = 0.3;
    double lnsqrts = std::log(sqrts);
    int pid;
    for (int j=0; j<4; j++){
        if (j<2) pid=4;
        else pid=5;
        std::ofstream f("Xsamples/"+S[j]+".dat");    
        for (int i=0; i<10000; i++){
            std::vector<fourvec> FS;
            std::vector<int> ID;
            X[j].sample({lnsqrts, T}, pid, FS, ID);
            if (FS.size()==3)
                f << ID[0] << " " << FS[0] << " "
                  << ID[1] << " " << FS[1] << " "
                  << ID[2] << " " << FS[2] << std::endl; 
        }
        f.close();
    }
}
