#include "Pythia8/Pythia.h"
#include <fstream>
#include <vector>
#include <string>
using namespace Pythia8;

int main(int argc, char* argv[]){
      
    std::ifstream fi(argv[1]);
    Pythia B2e;   
    B2e.readString("Random:setSeed = on");
    B2e.readString("Random:seed = 0");
    B2e.readString("ProcessLevel:all = off");
    B2e.readString("HadronLevel:all = on");
    B2e.readString("HadronLevel:decay = on");
    B2e.init();

    Pythia B2D;
    B2D.readString("Random:setSeed = on");
    B2D.readString("Random:seed = 0");
    B2D.readString("ProcessLevel:all = off");
    B2D.readString("HadronLevel:all = on");
    B2D.readString("HadronLevel:decay = on");
    B2D.readString("411:mayDecay = off");
    B2D.readString("421:mayDecay = off");
    B2D.readString("413:mayDecay = off");
    B2D.readString("423:mayDecay = off");
    B2D.init();


    Pythia B2JP;
    B2JP.readString("Random:setSeed = on");
    B2JP.readString("Random:seed = 0");
    B2JP.readString("ProcessLevel:all = off");
    B2JP.readString("HadronLevel:all = on");
    B2JP.readString("HadronLevel:decay = on");
    B2JP.readString("443:mayDecay = off");
    B2JP.init();




    // fill B hadrons
    std::string fname("e-"), f2(argv[1]);
    std::ofstream f(fname+f2);
    int line=0;
    std::string header;
    while (!fi.eof()){   
        if (line < 4) {
            std::getline(fi, header);
            f << header << std::endl;
            line++;
            continue;
        }
        line++;
        int j, pid;
        double px, py, pz, E, M, x, y, z, t, T, vx, vy, vz, px0, py0, pz0, w;
        fi >> j >> pid 
           >> px >> py >> pz >> E >> M
           >> x >> y   >> z  >> t 
           >> T >> vx  >> vy >> vz
           >> px0 >> py0 >> pz0 >> w;
        bool succeed = false;
        int apid = std::abs(pid); 

	 // output the original particle
         f  <<" "<< 4*j  <<" "<< pid
                   <<" "<< px << " "<< py  <<" "<< pz  <<" "<< E  <<" "<< M
                             <<" "<< x  <<" "<< y    <<" "<< z   <<" "<< t
                             <<" "<< T  <<" "<< vx   <<" "<< vy  <<" "<< vz
                             <<" "<< px0  <<" "<< py0  <<" "<< pz0  <<" "<< w<<" "<< std::endl;
        // pick B meson
        if (!((apid == 511)||(apid == 521)||(apid == 513)||(apid == 523))) 
            continue;
        // Force it into election, including both b->e and b->c->e
        while(!succeed){
            B2e.event.reset();
            B2e.event.append(pid, 85, 0, 0, px, py, pz, E, M);
            B2e.next();
	    //std::cout << "-----" << pythia.event.size() << std::endl;
            for (int i = 0; i < B2e.event.size(); ++i){
                auto p = B2e.event[i];
		//std::cout << p.id() << std::endl;
                if (p.idAbs()==11){
                    if (p.mother1()>0 && p.mother2()==0){
		      auto m1 = B2e.event[p.mother1()];
		      bool good = false;
		      if (m1.id()==pid) good = true;
		      else {
                         int abs = m1.idAbs();
                         if (abs==411||abs==421||abs==413||abs==423||
                             abs==431||abs==433||abs==4122) 
                         good=true;
		      } 
                      if ( good ) {
                          f  <<" "<< 4*j+1  <<" "<< p.id()
                             <<" "<< p.px() << " "<< p.py()  <<" "<< p.pz()  <<" "<< p.e()  <<" "<< p.m()
                             <<" "<< x  <<" "<< y    <<" "<< z   <<" "<< t 
                             <<" "<< T  <<" "<< vx   <<" "<< vy  <<" "<< vz
                             <<" "<< px0  <<" "<< py0  <<" "<< pz0  <<" "<< w<<" "<< std::endl;
                          succeed = true;
			  break;
                       }
		    }
                }
            }
        }
        // force it into D meson:
	succeed = false;
        while(!succeed){
            B2D.event.reset();
            B2D.event.append(pid, 85, 0, 0, px, py, pz, E, M);
            B2D.next();
            for (int i = 0; i < B2D.event.size(); ++i){
                auto p = B2D.event[i];
                bool isD = (p.idAbs()==411) || (p.idAbs()==421) 
                        || (p.idAbs()==413) || (p.idAbs()==423);
                if (isD){
                  succeed = true;
                  f  << 4*j+2  <<" "<< 540
                     <<" "<< p.px() << " " << p.py()  <<" "<< p.pz()  <<" "<< p.e()  <<" "<< p.m()
                     <<" "<< x  <<" "<< y    <<" "<< z   <<" "<< t
                     <<" "<< T  <<" "<< vx   <<" "<< vy  <<" "<< vz
                     <<" "<< px0  <<" "<< py0  <<" "<< pz0  <<" "<< w<<" "<< std::endl;
		  break;
                }
            }
        }

     	// force it into J/Psi:
        succeed = false;
        while(!succeed){
            B2JP.event.reset();
            B2JP.event.append(pid, 85, 0, 0, px, py, pz, E, M);
            B2JP.next();
            for (int i = 0; i < B2JP.event.size(); ++i){
                auto p = B2JP.event[i];
                bool isJP = p.idAbs()==443;
                if (isJP){
                  succeed = true;
                  f  << 4*j+3  <<" "<< 544
                     <<" "<< p.px() << " " << p.py()  <<" "<< p.pz()  <<" "<< p.e()  <<" "<< p.m()
                     <<" "<< x  <<" "<< y    <<" "<< z   <<" "<< t
                     <<" "<< T  <<" "<< vx   <<" "<< vy  <<" "<< vz
                     <<" "<< px0  <<" "<< py0  <<" "<< pz0  <<" "<< w<<" "<< std::endl;
		  break;
                }
            }
        }
    }
    return 0;

}

