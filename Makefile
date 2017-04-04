

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

tests/development.o: tests/development.cxx $(HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

bin/dev.app: tests/development.o
	@mkdir -p bin
	$(LD) $(LDFLAGS) ${ROOTGLIBS} ${ROOTLIBS} -o $@ $^

dev:bin/dev.app

tests/toymc.o: tests/toymc.cxx $(HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

bin/toymc.app: tests/toymc.o
	@mkdir -p bin
	$(LD) $(LDFLAGS) ${ROOTGLIBS} ${ROOTLIBS} -o $@ $^

toymc:bin/toymc.app

clean:
	@rm -f bin/*.app
	@rm -f tests/*.o


rootdict: MultiCumulants/LinkDef.h ToyMC/LinkDef.h
	rootcint -v4 -f MultiCumulants/cint_dictionary.cxx ToyMC/cint_dictionary.cxx -c $(HEADERS) $<
rootobj: MultiCumulants/cint_dictionary.cxx ToyMC/cint_dictionary.cxx
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(ROOTCFLAGS) $< -o MultiCumulants/cint_dictionary.o ToyMC/cint_dictionary.o

lib/MultiCumulants.so: MultiCumulants/cint_dictionary.o
	$(LD) $(SOFLAGS) $(LDFLAGS) $(ROOTLIBS) $^ -o $@

lib/ToyMC.so: ToyMC/cint_dictionary.o
	$(LD) $(SOFLAGS) $(LDFLAGS) $(ROOTLIBS) $^ -o $@

rootlib: lib/MultiCumulants.so lib/ToyMC.so
# DO NOT DELETE
