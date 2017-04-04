

NAME            := multicumulants
#http://semver.org/
VERSION         := 0.1.1


# ROOT
ROOT            := root
ROOTFLAGS       := -l -b -q 

ROOTCFLAGS    	= $(shell root-config --cflags)
ROOTLDFLAGS    	= $(shell root-config --ldflags)
ROOTLIBS      	= $(shell root-config --libs)
ROOTLDFLAGS    	= $(shell root-config --ldflags)
ROOT_ARCH      	= $(shell root-config --arch)


CXX             := g++ -c 
CXXFLAGS        := -std=c++14 -fPIC -Wall -Wextra -pedantic $(ROOTCFLAGS)
CPPFLAGS        := -I. -I./vendor/loguru

LD              := g++
LDFLAGS         := -Wl $(ROOTCFLAGS)
SOFLAGS         := -shared

#------------------------------------------------------------------------------
HEADERS 		:= MultiCumulants/QVector.h \
		           MultiCumulants/Subsets.h \
                           MultiCumulants/Algorithm.h \
                           ToyMC/ToyMCEvent.h \
                           ToyMC/ToyMCGenerator.h \
                           ToyMC/ToyMCParticle.h         

# Objects to compile into ToyMc binary
TOYMC 			:= tests/toymc.o

tests/development.o: tests/development.cxx $(HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

bin/dev.app: tests/development.o
	@mkdir -p bin
	$(LD) $(LDFLAGS) ${ROOTGLIBS} ${ROOTLIBS} -o $@ $^

.PHONY: dev
dev:bin/dev.app

# Produces .o files for ToyMC
$(TOYMC) : %.o: %.cxx
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

bin/toymc.app: $(TOYMC) ToyMC/cint_dictionary.o
	@mkdir -p bin
	$(LD) $(LDFLAGS) ${ROOTGLIBS} ${ROOTLIBS} -o $@ $^

.PHONY: toymc
toymc:bin/toymc.app

.PHONY: clean
clean:
	@rm -f bin/*.app
	@rm -f tests/*.o
	@rm -f ToyMC/*.o
	@rm -f ToyMC/cint_dictionary.*


ToyMC/cint_dictionary.cxx: ToyMC/LinkDef.h
	rootcint -v4 -f $@ -c ToyMC/ToyMCEvent.h ToyMC/ToyMCParticle.h ToyMC/ToyMCGenerator.h $<

ToyMC/cint_dictionary.o: ToyMC/cint_dictionary.cxx
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(ROOTCFLAGS) $< -o $@

# include the toymc.o because it has the implementation of the logger
lib/ToyMC.so: ToyMC/cint_dictionary.o $(TOYMC)
	$(LD) $(SOFLAGS) $(LDFLAGS) $(ROOTLIBS) $^ -o $@

.PHONY: rootlib
rootlib: lib/ToyMC.so
# DO NOT DELETE
