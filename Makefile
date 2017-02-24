

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
CXXFLAGS        := -std=c++14 -fPIC -Wall -Wextra -pedantic
CPPFLAGS        := -I.

LD              := g++
LDFLAGS         := -Wl
SOFLAGS         := -shared

#------------------------------------------------------------------------------
HEADERS 		:= MultiCumulants/QVector.h \
		           MultiCumulants/Subsets.h          

tests/development.o: tests/development.cxx $(HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

bin/dev.app: tests/development.o
	$(LD) $(LDFLAGS) -o $@ $^

dev:bin/dev.app


clean:
	@rm -f bin/dev.app
	@rm -f tests/*.o


rootdict: MultiCumulants/LinkDef.h
	rootcint -v4 -f MultiCumulants/cint_dictionary.cxx -c $(HEADERS) $<
rootobj: MultiCumulants/cint_dictionary.cxx
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(ROOTCFLAGS) $< -o MultiCumulants/cint_dictionary.o

lib/MultiCumulants.so: MultiCumulants/cint_dictionary.o
	$(LD) $(SOFLAGS) $(LDFLAGS) $(ROOTLIBS) $^ -o $@

rootlib: lib/MultiCumulants.so
