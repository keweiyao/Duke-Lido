# deal with pythia
PYTHIADIR   = /usr
PYTHIAFLAGS = -I$(PYTHIADIR)/include -O2 -pedantic -W -Wall -Wshadow -fPIC
PYTHIALIBS  =  -L$(PYTHIADIR)/lib -Wl,-rpath,$(PYTHIADIR) -lpythia8 -ldl

# deal with fastjet
FASTJETDIR  = /usr/lib
FASTJETARGS = `$(FASTJETDIR)/bin/fastjet-config --cxxflags --libs`

# LHAPDF if you want to replace Pythia's dummy default
LHAPDFARGS =
#LHAPDFARGS = -L$(HOME)/local/lib -lLHAPDF -I$(HOME)/local/include

#main72:
#	g++ -o main72 main72.cc $(PYTHIAFLAGS) $(PYTHIALIBS) $(FASTJETARGS) $(LHAPDFARGS) 
#main71:
#	g++ -o main71 main71.cc $(PYTHIAFLAGS) $(PYTHIALIBS) $(FASTJETARGS) $(LHAPDFARGS) 
#jet:
#	g++ -o Jet Jet_cs.cpp $(PYTHIAFLAGS) $(PYTHIALIBS) $(FASTJETARGS) $(LHAPDFARGS) 
jet_simu:
	g++ -std=c++11 -o jet_cross_section jet_cross_section.cpp $(PYTHIAFLAGS) $(PYTHIALIBS) $(FASTJETARGS) $(LHAPDFARGS) 
	g++ -std=c++11 -o jet_lund_diagram jet_lund_diagram.cpp $(PYTHIAFLAGS) $(PYTHIALIBS) $(FASTJETARGS) $(LHAPDFARGS) 
	g++ -std=c++11 -o jet_shape jet_shape.cpp $(PYTHIAFLAGS) $(PYTHIALIBS) $(FASTJETARGS) $(LHAPDFARGS) 
	g++ -std=c++11 -o jet_fragmentation jet_fragmentation.cpp $(PYTHIAFLAGS) $(PYTHIALIBS) $(FASTJETARGS) $(LHAPDFARGS) 
	g++ -std=c++11 -o dijet_asymmetry dijet_asymmetry.cpp $(PYTHIAFLAGS) $(PYTHIALIBS) $(FASTJETARGS) $(LHAPDFARGS) 
