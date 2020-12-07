#include <string>
#include <iostream>
#include <exception>

#include "Rate.h"
#include "simpleLogger.h"
#include "predefine.h"
#include "matrix_elements.h"
#include <fstream>
#include <vector>
#include <iterator>


void test_R22_q(void);
void test_R22_g(void);
void test_R22_hq(void);
void test_R23_q(void);
void test_R23_g(void);
void test_R23_hq(void);
void test_R12(void);

int main(){
    initialize_mD_and_scale(0, 2, 0.3, 4., 4.);
    LOG_INFO << "R22-q";
    test_R22_q();
    LOG_INFO << "R22-g";
    test_R22_g();
    LOG_INFO << "R22-hq";
    test_R22_hq();
    LOG_INFO << "R23-q";
    test_R23_q();
    LOG_INFO << "R23-g";
    test_R23_g();
    LOG_INFO << "R23-hq";
    test_R23_hq();
    LOG_INFO << "R12";
    test_R12();
}

void test_R22_q(void){
    Rate<HS2PP, 2, 2, double(*)(const double, void*)>
    r1("Boltzmann/qg2qg", "lido_setting.xml", dX_qg2qg),
    r2("Boltzmann/qq2qq", "lido_setting.xml", dX_qq2qq),
    r3("Boltzmann/qqbar2qqbar", "lido_setting.xml", dX_qqbar2qqbar_diff),
    r4("Boltzmann/qqbar2ccbar", "lido_setting.xml", dX_qqbar2qqbar_diff),
    r5("Boltzmann/qqbar2bbbar", "lido_setting.xml", dX_qqbar2qqbar_diff);
    std::vector<Rate<HS2PP, 2, 2, double(*)(const double, void*)> >
    R{r1, r2, r3, r4, r5};
    std::vector<std::string> S{"qg2qg", "qq2qq", "qqbar2qqbar",
                               "qqbar2ccbar", "qqbar2bbbar"};
    for (auto & r : R) {
     r.loadX("table.h5");
     r.load("table.h5");
    }
    double E= 20;
    double T = 0.3;
    double lnE = std::log(E);
    int pid=1;
    for (int j=0; j<5; j++){
        std::ofstream f("Rsamples/"+S[j]+".dat");    
        for (int i=0; i<10000; i++){
            std::vector<fourvec> FS;
            std::vector<int> ID;
            R[j].sample({lnE, T}, pid, FS, ID);
            if (FS.size()==2)
                f << ID[0] << " " << FS[0] << " "
                  << ID[1] << " " << FS[1] << std::endl; 
        }
        f.close();
    }
}

void test_R22_g(void){
    Rate<HS2PP, 2, 2, double(*)(const double, void*)>
    r1("Boltzmann/gg2gg", "lido_setting.xml", dX_gg2gg),
    r2("Boltzmann/gq2gq", "lido_setting.xml", dX_gq2gq),
    r3("Boltzmann/gg2qqbar", "lido_setting.xml", dX_gg2qqbar),
    r4("Boltzmann/gg2ccbar", "lido_setting.xml", dX_gg2qqbar),
    r5("Boltzmann/gg2bbbar", "lido_setting.xml", dX_gg2qqbar);
    std::vector<Rate<HS2PP, 2, 2, double(*)(const double, void*)> >
    R{r1, r2, r3, r4, r5};
    std::vector<std::string> S{"gg2gg", "gq2gq", "gg2qqbar",
                               "gg2ccbar", "gg2bbbar"};
    for (auto & r : R) {
     r.loadX("table.h5");
     r.load("table.h5");
    }
    double E= 20;
    double T = 0.3;
    double lnE = std::log(E);
    int pid=21;
    for (int j=0; j<5; j++){
        std::ofstream f("Rsamples/"+S[j]+".dat");    
        for (int i=0; i<10000; i++){
            std::vector<fourvec> FS;
            std::vector<int> ID;
            R[j].sample({lnE, T}, pid, FS, ID);
            if (FS.size()==2)
                f << ID[0] << " " << FS[0] << " "
                  << ID[1] << " " << FS[1] << std::endl; 
        }
        f.close();
    }
}

void test_R22_hq(void){
    Rate<HS2PP, 2, 2, double(*)(const double, void*)>
    r1("Boltzmann/cg2cg", "lido_setting.xml", dX_qg2qg),
    r2("Boltzmann/cq2cq", "lido_setting.xml", dX_qq2qq),
    r3("Boltzmann/bg2bg", "lido_setting.xml", dX_qg2qg),
    r4("Boltzmann/bq2bq", "lido_setting.xml", dX_qq2qq);
    std::vector<Rate<HS2PP, 2, 2, double(*)(const double, void*)> >
    R{r1, r2, r3, r4};
    std::vector<std::string> S{"cg2cg", "cq2cq", "bg2bg", "bq2bq"};
    for (auto & r : R) {
     r.loadX("table.h5");
     r.load("table.h5");
    }
    double E= 20;
    double T = 0.3;
    double lnE = std::log(E);
    int pid;
    for (int j=0; j<4; j++){
        pid=(j<2)? 4:5;
        std::ofstream f("Rsamples/"+S[j]+".dat");    
        for (int i=0; i<10000; i++){
            std::vector<fourvec> FS;
            std::vector<int> ID;
            R[j].sample({lnE, T}, pid, FS, ID);
            if (FS.size()==2)
                f << ID[0] << " " << FS[0] << " "
                  << ID[1] << " " << FS[1] << std::endl; 
        }
        f.close();
    }
}

void test_R23_q(void){
    Rate<HS2PPP, 2, 2, double(*)(const double *, void*)>
    r1("Boltzmann/qg2qgg", "lido_setting.xml", dX_qg2qgg),
    r2("Boltzmann/qq2qqg", "lido_setting.xml", dX_qq2qqg),
    r3("Boltzmann/qg2qqqbar", "lido_setting.xml", dX_qg2qqqbar),
    r4("Boltzmann/qg2qccbar", "lido_setting.xml", dX_qg2qqqbar),
    r5("Boltzmann/qg2qbbbar", "lido_setting.xml", dX_qg2qqqbar);
    std::vector<Rate<HS2PPP, 2, 2, double(*)(const double *, void*)> >
    R{r1, r2, r3, r4, r5};
    std::vector<std::string> S{"qg2qgg", "qq2qqg", "qg2qqqbar",
                               "qg2qccbar", "qg2qbbbar"};
    for (auto & r : R) {
     r.loadX("table.h5");
     r.load("table.h5");
    }
    double E= 100;
    double T = 0.3;
    double lnE = std::log(E);
    int pid=1;
    for (int j=0; j<5; j++){
        std::ofstream f("Rsamples/"+S[j]+".dat");    
        for (int i=0; i<10000; i++){
            std::vector<fourvec> FS;
            std::vector<int> ID;
            R[j].sample({lnE, T}, pid, FS, ID);
            if (FS.size()==3)
                f << ID[0] << " " << FS[0] << " "
                  << ID[1] << " " << FS[1] << " "
                  << ID[2] << " " << FS[2] << std::endl; 
        }
        f.close();
    }
}

void test_R23_g(void){
    Rate<HS2PPP, 2, 2, double(*)(const double *, void*)>
    r1("Boltzmann/gg2ggg", "lido_setting.xml", dX_gg2ggg),
    r2("Boltzmann/gq2gqg", "lido_setting.xml", dX_gq2gqg),
    r3("Boltzmann/gg2qgqbar", "lido_setting.xml", dX_gg2qgqbar),
    r4("Boltzmann/gg2cgcbar", "lido_setting.xml", dX_gg2qgqbar),
    r5("Boltzmann/gg2bgbbar", "lido_setting.xml", dX_gg2qgqbar),
    r6("Boltzmann/gq2qqqbar", "lido_setting.xml", dX_gq2qqqbar),
    r7("Boltzmann/gq2cqcbar", "lido_setting.xml", dX_gq2qqqbar),
    r8("Boltzmann/gq2bqbbar", "lido_setting.xml", dX_gq2qqqbar);
    std::vector<Rate<HS2PPP, 2, 2, double(*)(const double *, void*)> >
    R{r1, r2, r3, r4, r5, r6, r7, r8};
    std::vector<std::string> S{"gg2ggg", "gq2gqg", "gg2qgqbar", "gg2cgcbar",
                        "gg2bgbbar", "gq2qqqbar", "gq2cqcbarr", "gq2bqbbar"};
    for (auto & r : R) {
     r.loadX("table.h5");
     r.load("table.h5");
    }
    double E= 100;
    double T = 0.3;
    double lnE = std::log(E);
    int pid=21;
    for (int j=0; j<8; j++){
        std::ofstream f("Rsamples/"+S[j]+".dat");    
        for (int i=0; i<10000; i++){
            std::vector<fourvec> FS;
            std::vector<int> ID;
            R[j].sample({lnE, T}, pid, FS, ID);
            if (FS.size()==3)
                f << ID[0] << " " << FS[0] << " "
                  << ID[1] << " " << FS[1] << " "
                  << ID[2] << " " << FS[2] << std::endl; 
        }
        f.close();
    }
}


void test_R23_hq(void){
    Rate<HS2PPP, 2, 2, double(*)(const double *, void*)>
    r1("Boltzmann/cg2cgg", "lido_setting.xml", dX_qg2qgg),
    r2("Boltzmann/cq2cqg", "lido_setting.xml", dX_qq2qqg),
    r3("Boltzmann/bg2bgg", "lido_setting.xml", dX_qg2qgg),
    r4("Boltzmann/bq2bqg", "lido_setting.xml", dX_qq2qqg);
    std::vector<Rate<HS2PPP, 2, 2, double(*)(const double *, void*)> >
    R{r1, r2, r3, r4};
    std::vector<std::string> S{"cg2cgg", "cq2cqg", "bg2bgg", "bq2bqg"};
    for (auto & r : R) {
     r.loadX("table.h5");
     r.load("table.h5");
    }
    double E = 100;
    double T = 0.3;
    double lnE = std::log(E);
    int pid;
    for (int j=0; j<4; j++){
        pid = (j<2)? 4:5;
        std::ofstream f("Rsamples/"+S[j]+".dat");    
        for (int i=0; i<10000; i++){
            std::vector<fourvec> FS;
            std::vector<int> ID;
            R[j].sample({lnE, T}, pid, FS, ID);
            if (FS.size()==3)
                f << ID[0] << " " << FS[0] << " "
                  << ID[1] << " " << FS[1] << " "
                  << ID[2] << " " << FS[2] << std::endl; 
        }
        f.close();
    }
}


void test_R12(void){
    EffRate12<2, double(*)(const double*, void*)>
    r1("Boltzmann/q2qg", "lido_setting.xml", LGV_q2qg),
    r2("Boltzmann/c2cg", "lido_setting.xml", LGV_q2qg),
    r3("Boltzmann/b2bg", "lido_setting.xml", LGV_q2qg),
    r4("Boltzmann/g2gg", "lido_setting.xml", LGV_g2gg),
    r5("Boltzmann/g2qqbar", "lido_setting.xml", LGV_g2qqbar),
    r6("Boltzmann/g2ccbar", "lido_setting.xml", LGV_g2qqbar),
    r7("Boltzmann/g2bbbar", "lido_setting.xml", LGV_g2qqbar);

    std::vector<EffRate12<2, double(*)(const double*, void*)> >
    R{r1, r2, r3, r4, r5, r6, r7};
    std::vector<std::string> S{"q2qg", "c2cg", "b2bg",
                               "g2gg", "g2qqbar",  "g2ccbar",  "g2bbbar"};
    for (auto & r : R)  r.load("table.h5");
    double E= 100;
    double T = 0.3;
    double lnE = std::log(E);
    int pids[7] = {1,4,5,21,21,21,21};
    for (int j=0; j<7; j++){
        std::ofstream f("Rsamples/"+S[j]+".dat");    
        for (int i=0; i<10000; i++){
            std::vector<fourvec> FS;
            std::vector<int> ID;
            R[j].sample({lnE, T}, pids[j], FS, ID);
            if (FS.size()==2)
                f << ID[0] << " " << FS[0] << " "
                  << ID[1] << " " << FS[1] << std::endl; 
        }
        f.close();
    }
}
