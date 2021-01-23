#include "predefine.h"
#include <cmath>
#include <boost/math/tools/roots.hpp>
#include "simpleLogger.h"
#include "random.h"

char HS2PP[] = "HS2PP"; 
char HS2PPP[] = "HS2PPP"; 
char HSS2PP[] = "HSS2PP"; 

bool type1_warned = false;
bool type2_warned = false; 
bool type3_warned = false; 

const double c2pi = 2.*M_PI;
const double c4d9 = 4./9.;
const double c1d9 = 1./9.;
const double c8pi = 8.*M_PI;
const double c16pi = 16.*M_PI;
const double c48pi = 48.*M_PI;
const double c16pi2 = 16.*M_PI*M_PI;
const double c72pi2 = 72.*M_PI*M_PI;
const double c32pi3 = 32.*M_PI*M_PI*M_PI;
const double c64d9pi2 = 64./9.*M_PI*M_PI;
const double c256pi4 = 256.*std::pow(M_PI, 4);
const double fmc_to_GeV_m1 = 5.076;
int color_count = -1;
// number of color=3 (3*3-1 = 8 gluons), number of flavor=3, (u,d,s quark)
const int Nc = 3, nf = 3;
const double CF = 4./3.;
const double CA = 3.;
const double CF_over_CA = CF/CA;
const double TR = 0.5;
const double Mc = 1.3;
const double Mb = 4.2;

// the prefractor for gluon debye mass with Boltzmann statistics
// mD^2 = 8\pi*(Nc+nf)*alpha_s*T^2 ~ 15*alpha_s*T^2
// Note that if using quantum statistics, this will be
// mD^2 = 4\pi/3*(Nc+nf/2)*alpha_s*T^2 ~ 18*alpha_s*T^2
const double pf_g = 8./M_PI*(Nc + nf); // prefractor for gluon self energy^2

// For QCD coupling constant
// alpha_s = alpha_0 = 4\pi/(11Nc/3-2Nf/3)/log(Q^2/LambdaQCD^2)
const double alpha0 = 4.*M_PI/(11./3.*Nc - 2./3.*nf); // alpha_s(Q2 = e*Lambda2)
const double alpha_max = 1.0;
const double Lambda = 0.2; // [GeV] Lambda QCD = 0.2 GeV
const double Lambda2 = Lambda*Lambda; // [GeV^2] Lambda QCD squared
const double mu2_NP = Lambda2*std::exp(alpha0/alpha_max); // minimum cut on Q2, where alpha = alpha_0

double cut;
double scale;
double afix;
double Lido_Ecut;

Debye_mass * t_channel_mD2;
const double LPM_prefactor = 2.0;

int time_type;
bool Adiabatic_LPM;

void initialize_mD_and_scale(int _mD_type, double _mu, double _afix, double _theta, double _cut){
        cut = _cut;
        scale = _mu;
        afix = _afix;
    Lido_Ecut = _theta;
        t_channel_mD2 = new Debye_mass(_mD_type);
}

double pid2mass(int pid){
    int abspid = std::abs(pid);
    if (abspid==123||abspid==1||abspid==2||abspid==3||abspid==21){
        return 0.;
    }
    else if (abspid==4) return Mc;
    else if (abspid==5) return Mb;
    else return 0.;
}


void echo(void){
    std::cout << "mu\tafix\tQcut\ttheta\n";
    std::cout <<scale<<"\t"<<afix<<"\t"<<cut<<"\t"<<Lido_Ecut<<"\n";
}

Debye_mass::Debye_mass(const unsigned int _type):
TL(0.1), TH(1.0), NT(100), dT((TH-TL)/(NT-1.)),
type(_type), mD2(new double[NT])
{
        if (type==0) {
                LOG_INFO << "# leading order Debye mass";
                // type==0 use self-consistent Debye mass
                for (size_t i=0; i<NT; i++){
                        double T = TL+dT*i;
                        mD2[i] = pf_g*alpha_s(0., T)*T*T;
                }
        }
        if (type==1) {
                LOG_INFO << "# self-consistent Debye mass";
                // type==1 use self-consistent Debye mass
                for (size_t i=0; i<NT; i++){
                        double T = TL+dT*i;
                        size_t maxiter=200;
                        boost::math::tools::eps_tolerance<double> tol{
                         (std::numeric_limits<double>::digits * 3) / 4};
                        try{
                                auto result = boost::math::tools::toms748_solve(
                                        [&T](double x) {return pf_g*alpha_s(-x, T)*T*T - x;},
                                        0.001, 100., tol, maxiter);
                                mD2[i] = .5*(result.first + result.second);
                        }
                        catch (const std::domain_error&) {
                                throw std::domain_error{
                                "unable to calculate mD2"};
                        }
                }
        }
}

double Debye_mass::get_mD2(double T){
        if (T<TL) T=TL;
        if (T>=TH-dT) T=TH-dT;
        double x = (T-TL)/dT;
        size_t index = std::floor(x);
        double r = x-index;
        return (1.-r)*mD2[index] + r*mD2[index+1];
}


// splitting functions
double P_q2qg(double x){
    return CF*(1. + std::pow(1.-x, 2))/x;
}
double P_q2gq(double x){
    return CF*(1. + std::pow(x, 2))/(1-x);
}
double P_g2gg(double x){
    return CA*(1.+std::pow(x, 4)+std::pow(1.-x,4))/x/(1.-x);
}
double P_g2qqbar(double x){
    return ( x*x + std::pow(1.-x, 2) )*TR;
}

// split=1: q->q+g, colors = 1 - x + CF/CA * x^2
// split=2: g->g+g, colors = 1 - x + x^2
// split=3: g->q+qbar, colors = 1 - x*CA/CF + x^2*CA/CF
void formation_time(double & tauf, double & Q2,
           int idA, int idB, int idC,
           fourvec pA, fourvec pB, fourvec pC, 
           double T, fourvec p0){
    double mg2 = t_channel_mD2->get_mD2(T) / 2.;
    fourvec Ptot = p0*((pB.t()+pC.t())/p0.t());
    double Ptotabs = Ptot.pabs();
    fourvec nbar{1.,-Ptot.x()/Ptotabs,-Ptot.y()/Ptotabs,-Ptot.z()/Ptotabs};
    double Etot = Ptotabs + Ptot.t();
    double x = dot(pB, nbar)/Etot;
    if (x<0.||x>1.) { 
        Q2 = mg2; 
        LOG_WARNING << "x<0 or x>1 in formation time, x = "<<x ;
        LOG_INFO << p0;
        LOG_INFO << pA;
        LOG_INFO << pB;
        LOG_INFO << pC;
    }
    double M2 = x*(std::pow(pid2mass(idC), 2)+ mg2)
            + (1.-x)*(std::pow(pid2mass(idB), 2)+ mg2)
          - x*(1.-x)*(std::pow(pid2mass(idA), 2) + mg2);
    Q2 =  x*(measure_perp(p0, pC)).pabs2()
       + (1.-x)*(measure_perp(p0, pB)).pabs2()
       - x*(1.-x)*(measure_perp(p0, pA)).pabs2() 
       + M2;
    if (Q2<M2) { 
        Q2 = M2;
        //LOG_WARNING << "Q2<M2 in formationt time, set it by M2";
    }
    tauf = x*(1.-x)*Etot/Q2;
}

//=============running coupling========================
double alpha_s(double Q2, double T){
        if (afix > 0.) {
                return afix;
        }
        else{
                double screen_scale2 = std::pow(scale*M_PI*T, 2);
                double mu2;
                mu2 = std::max(std::abs(Q2), screen_scale2);

                if (mu2 <= mu2_NP) return 1.0;
                else return alpha0 / std::log(mu2/Lambda2);
        }
}

bool is_file_exist(std::string fileName)
{
    std::ifstream infile(fileName);
    return infile.good();
}

fourvec put_on_shell(fourvec p, int pid){
    double m = pid2mass(pid);
    p.a[0] = std::sqrt(p.pabs2()+m*m);
    return p;
}


particle make_parton(int pid, int col, int acol, 
                     double scale, 
                     fourvec pmu, 
                     particle & mother,
                     double T, std::vector<double> vcell){
    particle p;
    p.pid = pid;
    p.mass = pid2mass(pid);
    p.weight = mother.weight;
    p.p0 = pmu;
    p.p = pmu;
    p.x0 = mother.x;
    p.x = mother.x;
    p.mother_p = mother.p;
    p.T0 = T;
    p.Tf = T;
    p.is_virtual = mother.is_virtual;
    p.vcell.resize(3);
    p.vcell[0] = vcell[0];
    p.vcell[1] = vcell[1];
    p.vcell[2] = vcell[2];
    p.col = col;
    p.acol = acol;
    p.Q00 = scale;
    p.Q0 = scale;
    p.mfp0 = mother.mfp0;
    return p;
}

particle make_parton(int pid, int col, int acol, 
                     double scale, 
                     fourvec pmu, 
                     coordinate x
                     ){
    particle p;
    p.pid = pid;
    p.mass = pid2mass(pid);
    p.weight = 1;
    p.p0 = pmu;
    p.p = pmu;
    p.x0 = x;
    p.x = x;
    p.mother_p = pmu;
    p.T0 = 0;
    p.Tf = 0;
    p.is_virtual = false;
    p.vcell.resize(3);
    p.vcell[0] = 0;
    p.vcell[1] = 0;
    p.vcell[2] = 0;
    p.col = col;
    p.acol = acol;
    p.Q00 = scale;
    p.Q0 = scale;
    p.mfp0 = 0;
    return p;
}

bool is1to2(int prcid){ // 7
    return (prcid==14||prcid==15||prcid==16||prcid==17
      ||prcid==28||prcid==33||prcid==38);
}
bool is2to2(int prcid){ // 14
    return (prcid==1||prcid==2||prcid==3||prcid==4||prcid==5||
            prcid==18||prcid==19||prcid==20||prcid==21||prcid==22||
            prcid==29||prcid==30||prcid==34||prcid==35);
}
bool is2to2_nlo(int prcid){
    return (prcid==25||prcid==26||prcid==27||
            prcid==39||prcid==40||prcid==41);
}
bool is_classical(int prcid){
    return (prcid == 1 || prcid == 2 
                    || prcid == 18 || prcid == 19
                    || prcid == 29 || prcid == 30
                    || prcid == 34 || prcid == 35 );
}
bool is2to3(int prcid){ // 17
    return (prcid==6||prcid==7||prcid==8||prcid==9
           ||prcid==10||prcid==11||prcid==12||prcid==13
           ||prcid==23||prcid==24||prcid==25||prcid==26
           ||prcid==27||prcid==31||prcid==32||prcid==36
           ||prcid==37||prcid==39||prcid==40||prcid==41);
}
bool isPairProduction(int prcid){
   return (prcid==3||prcid==4||prcid==5
         ||prcid==8||prcid==9||prcid==10
         ||prcid==11||prcid==12||prcid==13
         ||prcid==15||prcid==16||prcid==17
         ||prcid==20||prcid==21||prcid==22
         ||prcid==25||prcid==26||prcid==27
         ||prcid==39||prcid==40||prcid==41);
}


int get_process_info(std::string p,
        std::vector<double> & _IM,
        std::vector<double> & _FM,
        std::vector<int> & _IT,
        std::vector<int> & _FT
){
    if (p=="gq2gq"){ // 1
        _IM.resize(2); _FM.resize(2); _IT.resize(2); _FT.resize(2);
        _IM[0]=0.; _IM[1]=0.; _FM[0]=0.; _FM[1]=0.;
        _IT[0]=21; _IT[1]=123; _FT[0]=21; _FT[1]=123;      
        return 1;  
    } 
    else if (p=="gg2gg"){ // 2
        _IM.resize(2); _FM.resize(2); _IT.resize(2); _FT.resize(2);
        _IM[0]=0.; _IM[1]=0.; _FM[0]=0.; _FM[1]=0.;
        _IT[0]=21; _IT[1]=21; _FT[0]=21; _FT[1]=21;   
        return 2;       
    } 
    else if (p=="gg2qqbar"){ // 3
        _IM.resize(2); _FM.resize(2); _IT.resize(2); _FT.resize(2);
        _IM[0]=0.; _IM[1]=0.; _FM[0]=0.; _FM[1]=0.;
        _IT[0]=21; _IT[1]=21; _FT[0]=123; _FT[1]=123;  
        return 3;        
    } 
    else if (p=="gg2ccbar"){ // 4
        _IM.resize(2); _FM.resize(2); _IT.resize(2); _FT.resize(2);
        _IM[0]=0.; _IM[1]=0.; _FM[0]=Mc; _FM[1]=Mc;
        _IT[0]=21; _IT[1]=21; _FT[0]=4; _FT[1]=4;   
        return 4;       
    } 
    else if (p=="gg2bbbar"){ // 5
        _IM.resize(2); _FM.resize(2); _IT.resize(2); _FT.resize(2);
        _IM[0]=0.; _IM[1]=0.; _FM[0]=Mb; _FM[1]=Mb;
        _IT[0]=21; _IT[1]=21; _FT[0]=5; _FT[1]=5;   
        return 5;       
    } 
    else if (p=="gq2gqg"){ // 6
        _IM.resize(2); _FM.resize(3); _IT.resize(2); _FT.resize(3);
        _IM[0]=0.; _IM[1]=0.; _FM[0]=0.; _FM[1]=0.; _FM[2]=0.;
        _IT[0]=21; _IT[1]=123; _FT[0]=21; _FT[1]=123; _FT[2]=21; 
        return 6;         
    } 
    else if (p=="gg2ggg"){ // 7
        _IM.resize(2); _FM.resize(3); _IT.resize(2); _FT.resize(3);
        _IM[0]=0.; _IM[1]=0.; _FM[0]=0.; _FM[1]=0.; _FM[2]=0.;
        _IT[0]=21; _IT[1]=21; _FT[0]=21; _FT[1]=21; _FT[2]=21;   
        return 7;       
    } 
    else if (p=="gq2qqqbar"){ // 8
        _IM.resize(2); _FM.resize(3); _IT.resize(2); _FT.resize(3);
        _IM[0]=0.; _IM[1]=0.; _FM[0]=0.; _FM[1]=0.; _FM[2]=0.;
        _IT[0]=21; _IT[1]=123; _FT[0]=123; _FT[1]=123; _FT[2]=123; 
        return 8;         
    } 
    else if (p=="gq2cqcbar"){ // 9
        _IM.resize(2); _FM.resize(3); _IT.resize(2); _FT.resize(3);
        _IM[0]=0.; _IM[1]=0.; _FM[0]=Mc; _FM[1]=0.; _FM[2]=Mc;
        _IT[0]=21; _IT[1]=123; _FT[0]=4; _FT[1]=123; _FT[2]=4;
        return 9;          
    } 
    else if (p=="gq2bqbbar"){ // 10
        _IM.resize(2); _FM.resize(3); _IT.resize(2); _FT.resize(3);
        _IM[0]=0.; _IM[1]=0.; _FM[0]=Mb; _FM[1]=0.; _FM[2]=Mb;
        _IT[0]=21; _IT[1]=123; _FT[0]=5; _FT[1]=123; _FT[2]=5;   
        return 10;       
    } 
    else if (p=="gg2qgqbar"){ // 11
        _IM.resize(2); _FM.resize(3); _IT.resize(2); _FT.resize(3);
        _IM[0]=0.; _IM[1]=0.; _FM[0]=0.; _FM[1]=0.; _FM[2]=0.;
        _IT[0]=21; _IT[1]=21; _FT[0]=123; _FT[1]=21; _FT[2]=123;   
        return 11;       
    } 
    else if (p=="gg2cgcbar"){ // 12
        _IM.resize(2); _FM.resize(3); _IT.resize(2); _FT.resize(3);
        _IM[0]=0.; _IM[1]=0.; _FM[0]=Mc; _FM[1]=0.; _FM[2]=Mc;
        _IT[0]=21; _IT[1]=21; _FT[0]=4; _FT[1]=21; _FT[2]=4;   
        return 12;       
    } 
    else if (p=="gg2bgbbar"){ // 13
        _IM.resize(2); _FM.resize(3); _IT.resize(2); _FT.resize(3);
        _IM[0]=0.; _IM[1]=0.; _FM[0]=Mb; _FM[1]=0.; _FM[2]=Mb;
        _IT[0]=21; _IT[1]=21; _FT[0]=5; _FT[1]=21; _FT[2]=5;   
        return 13;       
    } 
    else if (p=="g2gg"){ // 14
        _IM.resize(1); _FM.resize(2); _IT.resize(1); _FT.resize(2);
        _IM[0]=0.;  _FM[0]=0.; _FM[1]=0.;
        _IT[0]=21;  _FT[0]=21; _FT[1]=21;
        return 14;       
    } 
    else if (p=="g2qqbar"){ // 15
        _IM.resize(1); _FM.resize(2); _IT.resize(1); _FT.resize(2);
        _IM[0]=0.;  _FM[0]=0.; _FM[1]=0.;
        _IT[0]=21;  _FT[0]=123; _FT[1]=123;
        return 15;       
    } 
    else if (p=="g2ccbar"){ // 16
        _IM.resize(1); _FM.resize(2); _IT.resize(1); _FT.resize(2);
        _IM[0]=0.;  _FM[0]=Mc; _FM[1]=Mc;
        _IT[0]=21;  _FT[0]=4; _FT[1]=4;
        return 16;       
    } 
    else if (p=="g2bbbar"){ // 17
        _IM.resize(1); _FM.resize(2); _IT.resize(1); _FT.resize(2);
        _IM[0]=0.;  _FM[0]=Mb; _FM[1]=Mb;
        _IT[0]=21;  _FT[0]=5; _FT[1]=5;
        return 17;       
    } 
    else if (p=="qq2qq"){ // 18
        _IM.resize(2); _FM.resize(2); _IT.resize(2); _FT.resize(2);
        _IM[0]=0.; _IM[1]=0.; _FM[0]=0.; _FM[1]=0.;
        _IT[0]=123; _IT[1]=123; _FT[0]=123; _FT[1]=123;
        return 18;       
    } 
    else if (p=="qg2qg"){ // 19
        _IM.resize(2); _FM.resize(2); _IT.resize(2); _FT.resize(2);
        _IM[0]=0.; _IM[1]=0.; _FM[0]=0.; _FM[1]=0.;
        _IT[0]=123; _IT[1]=21; _FT[0]=123; _FT[1]=21;
        return 19;       
    } 
    else if (p=="qqbar2qqbar"){ // 20
        _IM.resize(2); _FM.resize(2); _IT.resize(2); _FT.resize(2);
        _IM[0]=0.; _IM[1]=0.; _FM[0]=0.; _FM[1]=0.;
        _IT[0]=123; _IT[1]=123; _FT[0]=123; _FT[1]=123;
        return 20;       
    } 
    else if (p=="qqbar2ccbar"){ // 21
        _IM.resize(2); _FM.resize(2); _IT.resize(2); _FT.resize(2);
        _IM[0]=0.; _IM[1]=0.; _FM[0]=Mc; _FM[1]=Mc;
        _IT[0]=123; _IT[1]=123; _FT[0]=4; _FT[1]=4;
        return 21;       
    } 
    else if (p=="qqbar2bbbar"){ // 22
        _IM.resize(2); _FM.resize(2); _IT.resize(2); _FT.resize(2);
        _IM[0]=0.; _IM[1]=0.; _FM[0]=Mb; _FM[1]=Mb;
        _IT[0]=123; _IT[1]=123; _FT[0]=5; _FT[1]=5;
        return 22;       
    } 
    else if (p=="qq2qqg"){ // 23
        _IM.resize(2); _FM.resize(3); _IT.resize(2); _FT.resize(3);
        _IM[0]=0.; _IM[1]=0.; _FM[0]=0.; _FM[1]=0.; _FM[2]=0.;
        _IT[0]=123; _IT[1]=123; _FT[0]=123; _FT[1]=123; _FT[2]=21;
        return 23;       
    } 
    else if (p=="qg2qgg"){ // 24
        _IM.resize(2); _FM.resize(3); _IT.resize(2); _FT.resize(3);
        _IM[0]=0.; _IM[1]=0.; _FM[0]=0.; _FM[1]=0.; _FM[2]=0.;
        _IT[0]=123; _IT[1]=21; _FT[0]=123; _FT[1]=21; _FT[2]=21;
        return 24;       
    } 
    else if (p=="qg2qqqbar"){ // 25
        _IM.resize(2); _FM.resize(3); _IT.resize(2); _FT.resize(3);
        _IM[0]=0.; _IM[1]=0.; _FM[0]=0.; _FM[1]=0.; _FM[2]=0.;
        _IT[0]=123; _IT[1]=21; _FT[0]=123; _FT[1]=123; _FT[2]=123;
        return 25;       
    } 
    else if (p=="qg2qccbar"){ // 26
        _IM.resize(2); _FM.resize(3); _IT.resize(2); _FT.resize(3);
        _IM[0]=0.; _IM[1]=0.; _FM[0]=0.; _FM[1]=Mc; _FM[2]=Mc;
        _IT[0]=123; _IT[1]=21; _FT[0]=123; _FT[1]=4; _FT[2]=4;
        return 26;       
    } 
    else if (p=="qg2qbbbar"){ // 27
        _IM.resize(2); _FM.resize(3); _IT.resize(2); _FT.resize(3);
        _IM[0]=0.; _IM[1]=0.; _FM[0]=0.; _FM[1]=Mb; _FM[2]=Mb;
        _IT[0]=123; _IT[1]=21; _FT[0]=123; _FT[1]=5; _FT[2]=5;
        return 27;       
    } 
    else if (p=="q2qg"){ // 28
        _IM.resize(1); _FM.resize(2); _IT.resize(1); _FT.resize(2);
        _IM[0]=0.;_FM[0]=0; _FM[1]=0.;
        _IT[0]=123; _FT[0]=123; _FT[1]=21;
        return 28;       
    } 
    else if (p=="cq2cq"){ // 29
        _IM.resize(2); _FM.resize(2); _IT.resize(2); _FT.resize(2);
        _IM[0]=Mc; _IM[1]=0.; _FM[0]=Mc; _FM[1]=0.;
        _IT[0]=4; _IT[1]=123; _FT[0]=4; _FT[1]=123;
        return 29;       
    } 
    else if (p=="cg2cg"){ // 30
        _IM.resize(2); _FM.resize(2); _IT.resize(2); _FT.resize(2);
        _IM[0]=Mc; _IM[1]=0.; _FM[0]=Mc; _FM[1]=0.;
        _IT[0]=4; _IT[1]=21; _FT[0]=4; _FT[1]=21;
        return 30;       
    } 
    else if (p=="cq2cqg"){ // 31
        _IM.resize(2); _FM.resize(3); _IT.resize(2); _FT.resize(3);
        _IM[0]=Mc; _IM[1]=0.; _FM[0]=Mc; _FM[1]=0.; _FM[2]=0.;
        _IT[0]=4; _IT[1]=123; _FT[0]=4; _FT[1]=123; _FT[2]=21;
        return 31;             
    } 
    else if (p=="cg2cgg"){ // 32
        _IM.resize(2); _FM.resize(3); _IT.resize(2); _FT.resize(3);
        _IM[0]=Mc; _IM[1]=0.; _FM[0]=Mc; _FM[1]=0.; _FM[2]=0.;
        _IT[0]=4; _IT[1]=21; _FT[0]=4; _FT[1]=21; _FT[2]=21;
        return 32;             
    } 
    else if (p=="c2cg"){ // 33
        _IM.resize(1); _FM.resize(2); _IT.resize(1); _FT.resize(2);
        _IM[0]=Mc;  _FM[0]=Mc; _FM[1]=0.; 
        _IT[0]=4; _FT[0]=4; _FT[1]=21; 
        return 33;             
    } 
    else if (p=="bq2bq"){ // 34
        _IM.resize(2); _FM.resize(2); _IT.resize(2); _FT.resize(2);
        _IM[0]=Mb; _IM[1]=0.; _FM[0]=Mb; _FM[1]=0.;
        _IT[0]=5; _IT[1]=123; _FT[0]=5; _FT[1]=123;
        return 34;       
    } 
    else if (p=="bg2bg"){ // 35
        _IM.resize(2); _FM.resize(2); _IT.resize(2); _FT.resize(2);
        _IM[0]=Mb; _IM[1]=0.; _FM[0]=Mb; _FM[1]=0.;
        _IT[0]=5; _IT[1]=21; _FT[0]=5; _FT[1]=21;
        return 35;       
    } 
    else if (p=="bq2bqg"){ // 36
        _IM.resize(2); _FM.resize(3); _IT.resize(2); _FT.resize(3);
        _IM[0]=Mb; _IM[1]=0.; _FM[0]=Mb; _FM[1]=0.; _FM[2]=0.;
        _IT[0]=5; _IT[1]=123; _FT[0]=5; _FT[1]=123; _FT[2]=21;
        return 36;             
    } 
    else if (p=="bg2bgg"){ // 37
        _IM.resize(2); _FM.resize(3); _IT.resize(2); _FT.resize(3);
        _IM[0]=Mb; _IM[1]=0.; _FM[0]=Mb; _FM[1]=0.; _FM[2]=0.;
        _IT[0]=5; _IT[1]=21; _FT[0]=5; _FT[1]=21; _FT[2]=21;
        return 37;             
    } 
    else if (p=="b2bg"){ // 38
        _IM.resize(1); _FM.resize(2); _IT.resize(1); _FT.resize(2);
        _IM[0]=Mb;  _FM[0]=Mb; _FM[1]=0.; 
        _IT[0]=5; _FT[0]=5; _FT[1]=21; 
        return 38;             
    } 
    else if (p=="gg2gqqbar"){ // 39
        _IM.resize(2); _FM.resize(3); _IT.resize(2); _FT.resize(3);
        _IM[0]=0.; _IM[1]=0.; _FM[0]=0.; _FM[1]=0.; _FM[2]=0.;
        _IT[0]=21; _IT[1]=21; _FT[0]=21; _FT[1]=123; _FT[2]=123;
        return 39;       
    } 
    else if (p=="gg2gccbar"){ // 40
        _IM.resize(2); _FM.resize(3); _IT.resize(2); _FT.resize(3);
        _IM[0]=0.; _IM[1]=0.; _FM[0]=0.; _FM[1]=Mc; _FM[2]=Mc;
        _IT[0]=21; _IT[1]=21; _FT[0]=21; _FT[1]=4; _FT[2]=4;
        return 40;       
    } 
    else if (p=="gg2gbbbar"){ // 41
        _IM.resize(2); _FM.resize(3); _IT.resize(2); _FT.resize(3);
        _IM[0]=0.; _IM[1]=0.; _FM[0]=0.; _FM[1]=Mb; _FM[2]=Mb;
        _IT[0]=21; _IT[1]=21; _FT[0]=21; _FT[1]=5; _FT[2]=5;
        return 41;       
    } 
}

void assign_1to2_pid(int process_id, 
                     int incoming_hard_id, 
                     std::vector<int> & pids){ // 7 channels
    pids.resize(2);
    switch (process_id){
        case 14: //g2gg
            pids[0] = 21; 
            pids[1] = 21;
            break;
        case 15: //g2qqbar
            pids[0] = Srandom::sample_flavor(nf); 
            pids[1] = -pids[0];
            break;
        case 16: //g2ccbar
            pids[0] = Srandom::binary_choice()? 4:-4;
            pids[1] = -pids[0];
            break;
        case 17: //g2bbbar
            pids[0] = Srandom::binary_choice()? 5:-5;
            pids[1] = -pids[0];
            break;
        case 28: //q2qg
            pids[0] = incoming_hard_id;
            pids[1] = 21;
            break;
        case 33: //c2cg
            pids[0] = incoming_hard_id;
            pids[1] = 21;
            break;
        case 38: //b2bg
            pids[0] = incoming_hard_id;
            pids[1] = 21;
            break;
        default:
            LOG_FATAL << "Assign pid: process_id " << process_id << " is not a known 1-2-body process";
            exit(-1);
            break;
    }
}

void assign_2to2_pid(int process_id, 
                     int incoming_hard_id, 
                     std::vector<int>& pids){  // 14 channels
    pids.resize(2);
    switch (process_id){
        case 1: // gq2gq
            pids[0] = incoming_hard_id; 
            pids[1] = Srandom::sample_flavor(nf);
            break;
        case 2: // gg2gg
            pids[0] = incoming_hard_id; 
            pids[1] = 21;
            break;
        case 3: // gg2qqbar
            pids[0] = Srandom::sample_flavor(nf);
            pids[1] = -pids[0];
            break;
        case 4: // gg2ccbar
            pids[0] = Srandom::binary_choice()? 4:-4;
            pids[1] = -pids[0];
            break;
        case 5: // gg2bbbar
            pids[0] = Srandom::binary_choice()? 5:-5;
            pids[1] = -pids[0];
            break;
        case 18: // qq2qq
            pids[0] = incoming_hard_id;
            pids[1] = Srandom::sample_flavor(nf);
            break;   
        case 19: // qg2qg
            pids[0] = incoming_hard_id;
            pids[1] = 21;
            break; 
        case 20: // qqbar2qqbar, diff
            int newpid;
            do{ 
                newpid = Srandom::sample_flavor(nf);
            }while(std::abs(incoming_hard_id)==std::abs(newpid));
            pids[0] = newpid;
            pids[1] = -newpid;
            break;   
        case 21: // qqbar2ccbar
            pids[0] = Srandom::binary_choice()? 4:-4;
            pids[1] = -pids[0];
            break;     
        case 22: // qqbar2bbbar
            pids[0] = Srandom::binary_choice()? 5:-5;
            pids[1] = -pids[0];
            break;    
        case 29: // cq2cq       
            pids[0] = incoming_hard_id;
            pids[1] = Srandom::sample_flavor(nf);
            break;     
        case 30: // cg2cg       
            pids[0] = incoming_hard_id;
            pids[1] = 21;
            break;   
        case 34: // bq2bq       
            pids[0] = incoming_hard_id;
            pids[1] = Srandom::sample_flavor(nf);
            break;     
        case 35: // bg2bg       
            pids[0] = incoming_hard_id;
            pids[1] = 21;
            break;     
        default:
            LOG_FATAL << "Assign pid: process_id " << process_id << " is not a known 2-2-body process";
            exit(-1);
            break;
        }
}

void assign_2to3_pid(int process_id,
                     int incoming_hard_id, 
                     std::vector<int> & pids, 
                     bool forward){
    pids.resize(3);  // 17 channels
    switch (process_id){
        case 6: // gq2gqg
            pids[0] = incoming_hard_id; 
            pids[1] = Srandom::sample_flavor(nf);
            pids[2] = 21;
            break;
        case 7: // gg2ggg
            pids[0] = incoming_hard_id; 
            pids[1] = 21;
            pids[2] = 21;
            break;
        case 8: // gq2qqqbar
            if (forward){
                pids[0] = Srandom::sample_flavor(nf); 
                pids[1] = Srandom::sample_flavor(nf);
                pids[2] = -pids[0];
            }else LOG_FATAL<< "process 8 \"g(F)+q(B)->q+q+qbar\" cannot be backward";
            break;  
        case 9: // gq2cqcbar
            if (forward){
                pids[0] = Srandom::binary_choice()? 4:-4;
                pids[1] = Srandom::sample_flavor(nf);
                pids[2] = -pids[0];
            }else LOG_FATAL<< "process 9 \"g(F)+q(B)->c+q+cbar\" cannot be backward";
            break;  
        case 10: // gq2bqbbar
            if (forward){
                pids[0] = Srandom::binary_choice()? 5:-5;
                pids[1] = Srandom::sample_flavor(nf);
                pids[2] = -pids[0];
            }else LOG_FATAL<< "process 10 \"g(F)+q(B)->b+q+bbar\" cannot be backward";
            break;  
        case 11: // gg2qgqbar
            if (forward){
                pids[0] = Srandom::sample_flavor(nf);
                pids[1] = 21;
                pids[2] = -pids[0];
            }else LOG_FATAL<< "process 11 \"g(F)+g(B)->q+g+qbar\" cannot be backward";
            break;  
        case 12: // gg2cgcbar
            if (forward){
                pids[0] = Srandom::binary_choice()? 4:-4;
                pids[1] = 21;
                pids[2] = -pids[0];
            }else LOG_FATAL<< "process 12 \"g(F)+g(B)->c+g+cbar\" cannot be backward";
            break;  
        case 13: // gg2bgbbar
            if (forward){
                pids[0] = Srandom::binary_choice()? 5:-5;
                pids[1] = 21;
                pids[2] = -pids[0];
            }else LOG_FATAL<< "process 13 \"g(F)+g(B)->b+g+bbar\" cannot be backward";
            break;  
        case 23: // qq2qqg 
            pids[0] = incoming_hard_id;
            pids[1] = Srandom::sample_flavor(nf);
            pids[2] = 21;
            break; 
        case 24: // qg2qgg
            pids[0] = incoming_hard_id;
            pids[1] = 21;
            pids[2] = 21;
            break; 
        case 25: // qg2qqqbar
            if (forward) LOG_FATAL<< "process 25 \"q(F)+g(B)->q+q+qbar\" cannot be forward";
            else{
                pids[0] = incoming_hard_id;
                pids[1] = Srandom::sample_flavor(nf);
                pids[2] = -pids[1];
            }
            break; 
        case 26: // qg2qccbar
            if (forward) LOG_FATAL<< "process 26 \"q(F)+g(B)->q+c+cbar\" cannot be forward";
            else{
                pids[0] = incoming_hard_id;
                pids[1] = Srandom::binary_choice()? 4:-4;
                pids[2] = -pids[1];
            }
            break; 
        case 27: // qg2qbbbar
            if (forward) LOG_FATAL<< "process 27 \"q(F)+g(B)->q+b+bbar\" cannot be forward";
            else{
                pids[0] = incoming_hard_id;
                pids[1] = Srandom::binary_choice()? 5:-5;
                pids[2] = -pids[1];
            }
            break; 
        case 31: // cq2cqg
            pids[0] = incoming_hard_id;
            pids[1] = Srandom::sample_flavor(nf);
            pids[2] = 21;
            break; 
        case 32: // cg2cgg
            pids[0] = incoming_hard_id;
            pids[1] = 21;
            pids[2] = 21;
            break; 
        case 36: // bq2bqg
            pids[0] = incoming_hard_id;
            pids[1] = Srandom::sample_flavor(nf);
            pids[2] = 21;
            break; 
        case 37: // bg2bgg
            pids[0] = incoming_hard_id;
            pids[1] = 21;
            pids[2] = 21;
            break; 
        case 39: // gg2gqqbar
            if (forward) LOG_FATAL<< "process 39 \"q(F)+g(B)->q+q+qbar\" cannot be forward";
            else{
                pids[0] = incoming_hard_id;
                pids[1] = Srandom::sample_flavor(nf);
                pids[2] = -pids[1];
            }
            break; 
        case 40: // gg2gccbar
            if (forward) LOG_FATAL<< "process 40 \"q(F)+g(B)->q+c+cbar\" cannot be forward";
            else{
                pids[0] = incoming_hard_id;
                pids[1] = Srandom::binary_choice()? 4:-4;
                pids[2] = -pids[1];
            }
            break; 
        case 41: // gg2gbbbar
            if (forward) LOG_FATAL<< "process 41 \"q(F)+g(B)->q+b+bbar\" cannot be forward";
            else{
                pids[0] = incoming_hard_id;
                pids[1] = Srandom::binary_choice()? 5:-5;
                pids[2] = -pids[1];
            }
            break; 
        default:
            LOG_FATAL << "Assign pid: process_id " << process_id << " is not a known 2-3-body process";
            exit(-1);
            break;
        }
}


void assign_2to2_color(int process_id, 
                       int pidA, int pid1, int pid2,
                       int Ca, int aCa,
                       int & C1, int & aC1,
                       int & C2, int & aC2
                       ){
    if (process_id==1){ // gq2gq
        if (pid2>0){
            C1 = Ca; aC1 = color_count; color_count ++;
            C2 = aC1; aC2 = 0;

        }else{
            C1 = color_count; aC1 = aCa; color_count ++;
            C2 = 0; aC2 = C1;
        }
    }
    else if (process_id==2){ // gg2gg
        if (Srandom::binary_choice()){
            C1 = Ca; aC1 = color_count; color_count ++;
            C2 = aC1; aC2 = color_count; color_count ++;
        }else{
            C1 = color_count; aC1 = aCa; color_count ++;
            C2 = color_count; aC2 = C1; color_count ++;
        }
    }
    else if (process_id==3 || process_id==4
             || process_id==5 ){ // gg2qqbar, ccbar, bbbar
        if (pid1>0){
            if (Srandom::binary_choice()){   
                C1 = Ca; aC1 = 0;
                C2 = 0; aC2 = color_count; color_count ++;
            }else{
                C1 = color_count; aC1 = 0; color_count ++;
                C2 = 0; aC2 = aCa; 
            }
        }
        else{
            if (Srandom::binary_choice()){   
                C1 = 0; aC1 = aCa;
                C2 = color_count; aC2 = 0; color_count ++;
            }else{
                C1 = 0; aC1 = color_count; color_count ++;
                C2 = Ca; aC2 = 0; 
            }
        }
    }
    else if (process_id == 18 || process_id == 29
             || process_id == 34){ // qq2qq
        if (pid1>0 && pid2>0){  
            C1 = color_count; aC1 = 0; color_count ++;
            C2 = Ca; aC2 = 0;
        }
        if (pid1<0 && pid2<0){  
            C1 = 0; aC1 = color_count; color_count ++;
            C2 = 0; aC2 = aCa;
        }
        if (pid1>0 && pid2<0){  
            C1 = color_count; aC1 = 0; color_count ++;
            C2 = 0; aC2 = C1;
        }
        if (pid1<0 && pid2>0){  
            C1 = 0; aC1 = color_count; color_count ++;
            C2 = aC1; aC2 = 0;
        }
    }
    else if (process_id == 19 || process_id == 30
             || process_id == 35) { // qg2qg
        if (pid1>0){
            C1 = color_count; aC1 = 0; color_count++;
            C2 = color_count; aC2 = C1; color_count++;
        }
        if (pid1<0){
            C1 = 0; aC1 = color_count; color_count++;
            C2 = aC1; aC2 = color_count; color_count++;
        }
    }
    else if (process_id==20 || process_id==21 
             || process_id==22){  // qq2qqbar, ccbar, bbbar
        if (pidA>0 && pid1>0){
            C1 = Ca; aC1 = 0;
            C2 = 0; aC2 = color_count; color_count++;
        }
        if (pidA<0 && pid1<0){
            C1 = 0; aC1 = aCa;
            C2 = color_count; aC2 = 0; color_count++;
        }   
        if (pidA>0 && pid1<0){
            C1 = 0; aC1 = color_count; color_count++;
            C2 = Ca; aC2 = 0; 
        }    
        if (pidA<0 && pid1>0){
            C1 = color_count; aC1 = 0; color_count++;
            C2 = 0; aC2 = aCa; 
        }    
    }else{
        LOG_FATAL << "Wrong process id when assiging color";
        exit(-1);
    }
}


void assign_n2np1_color(
                       int pidA, int pid1, int pid2,
                       int Ca, int aCa,
                       int & C1, int & aC1,
                       int & C2, int & aC2
                       ){  
    if (pidA==21 && pid1==21 && pid2==21){ // g->g+g
        if (Srandom::binary_choice()){
            C1 = Ca; aC1 = color_count; color_count ++;
            C2 = aC1; aC2 = aCa;
        }else{
            C1 = color_count; aC1 = aCa; color_count ++;
            C2 = Ca; aC2 = C1;
        }
    }
    else if (pidA==21&&pid1!=21&&pid2!=21){ //gH->q+qbar
        if (pid1>0){
            C1 = Ca; aC1 = 0;
            C2 = 0; aC2 = aCa;
        }
        else{
            C1 = 0; aC1 = aCa;
            C2 = Ca; aC2 = 0;
        }
    }
    else if (pidA!=21&&pid1!=21&&pid2==21
             ){ // q -> q+g
        if (pid1>0){
            C1 = color_count; aC1 = 0; color_count ++;
            C2 = Ca; aC2 = C1;
        }
        else{
            C1 = 0; aC1 = color_count; color_count ++;
            C2 = aC1; aC2 = aCa;
        }
    }else{
        LOG_FATAL << "Wrong process id when assiging color in n->n+1";
        LOG_FATAL << pidA << " " << pid1 << " " << pid2;
        exit(-1);
    }
}



void assign_2to2_nlo_color(int process_id, 
                       int pidA, int pid1, int pid2, int pid3,
                       int Ca, int aCa,
                       int & C1, int & aC1,
                       int & C2, int & aC2,
                       int & C3, int & aC3
                       ){
    if (process_id==25 || process_id==26 || process_id==27){ // qg2qqqbar
        if (pidA>0){
            C1 = color_count; aC1 = 0; color_count ++;
            if (pid2<0){
                C2 = 0; aC2 = C1;
                C3 = color_count, aC3 = 0; color_count ++;
            }else{
                C3 = 0; aC3 = C1;
                C2 = color_count, aC2 = 0; color_count ++;
            }
        }else{
            C1 = 0; aC1 = color_count; color_count ++;
            if (pid2>0){
                C2 = aC1; aC2 = 0;
                C3 = 0, aC3 = color_count; color_count ++;
            }else{
                C3 = aC1; aC3 = 0;
                C2 = 0, aC2 = color_count; color_count ++;
            }
        }
    }
    else if (process_id==39 || process_id==40 || process_id==41){ // gg2gqqbar
        if (Srandom::binary_choice()){
            C1 = Ca; aC1 = color_count; color_count ++;
            if (pid2>0){
                C2 = aC1; aC2 = 0;
                C3 = 0, aC3 = color_count; color_count ++;
            }else{
                C3 = aC1; aC3 = 0;
                C2 = 0, aC2 = color_count; color_count ++;
            }
        }else{
            C1 = color_count; aC1 = aCa; color_count ++;
            if (pid2<0){
                C2 = 0; aC2 = C1;
                C3 = color_count, aC3 = 0; color_count ++;
            }else{
                C3 = 0; aC3 = C1;
                C2 = color_count, aC2 = 0; color_count ++;
            }
        }
    }
    else{
        LOG_FATAL << "Wrong process id when assiging color";
        exit(-1);
    }
}


void output_oscar(const std::vector<particle> plist,
                  int abspid, std::string fname){
    // output OSCAR Format
    int Nparticles = 0;
    for (auto &p : plist)
    {
        if (std::abs(p.pid) == abspid)
            Nparticles++;
    }

    std::ofstream f(fname);
    f << "OSC1997A\n";
    f << "final_id_p_x\n";
    f << "    lbt  1.0alpha   208    82   208    82   aacm  0.1380E+04        1\n";
    f << "        1  " << std::setw(10)
      << Nparticles << "    0.001    0.001    1    1       1\n";

    int i = 0;
    for (auto &p : plist){
        if (std::abs(p.pid) == abspid){
            f << std::setw(10) << std::setfill(' ') << i << "  "      // particle index, i10,2x
              << std::setw(10) << std::setfill(' ') << p.pid << "  "; // particle id, i10,2x
            f << ff(p.p.x()) << "  "
              << ff(p.p.y()) << "  "
              << ff(p.p.z()) << "  "
              << ff(p.p.t()) << "  "
              << ff(p.mass) << "  "
              << ff(p.x.x1() / fmc_to_GeV_m1) << "  "
              << ff(p.x.x2() / fmc_to_GeV_m1) << "  "
              << ff(p.x.x3() / fmc_to_GeV_m1) << "  "
              << ff(p.x.x0() / fmc_to_GeV_m1) << "  "
              << ff(p.Tf) << "  "
              << ff(p.vcell[0]) << "  "
              << ff(p.vcell[1]) << "  "
              << ff(p.vcell[2]) << "  "
              << ff(p.p0.x()) << "  "
              << ff(p.p0.y()) << "  "
              << ff(p.p0.z()) << "  "
              << ff(p.weight) << "\n";
        }
        i++;
    }
}


