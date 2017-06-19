

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


CXX             := g++ -c -g 
CXXFLAGS        := -std=c++14 -fPIC -Wall -Wextra -pedantic $(ROOTCFLAGS)
CPPFLAGS        := -I. -I./vendor/loguru

LD              := g++
LDFLAGS         := -Wl $(ROOTCFLAGS)
SOFLAGS         := -shared

#------------------------------------------------------------------------------
HEADERS 		:= MultiCumulants/QVector.h \
		           MultiCumulants/Subsets.h \
                           MultiCumulants/Algorithm.h \
		           MultiCumulants/Correlator.h \
                           ToyMC/ToyMCEvent.h \
                           ToyMC/ToyMCGenerator.h \
                           ToyMC/ToyMCParticle.h \
                           correlations/Correlator.hh \
                           correlations/FromQVector.hh \
                           correlations/NestedLoops.hh \
                           correlations/QVector.hh \
                           correlations/Result.hh \
                           correlations/Types.hh \
                           correlations/closed/FromQVector.hh \
                           correlations/recurrence/FromQVector.hh \
                           correlations/recursive/FromQVector.hh \
                           correlations/recursive/NestedLoops.hh


# Objects to compile into ToyMc binary
TOYMC 			:= tests/toymc.o

# this works for make v 3.8 and newer
.DEFAULT_GOAL := toymc

###########################################################################
# TOY MonteCarlo binary
# Produces .o files for ToyMC
$(TOYMC) : %.o: %.cxx
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

bin/toymc.app: $(TOYMC) ToyMC/cint_dictionary.o
	@mkdir -p bin
	$(LD) $(LDFLAGS) ${ROOTGLIBS} ${ROOTLIBS} -o $@ $^

.PHONY: toymc
toymc:bin/toymc.app

toymc_depend: .depend_toymc

.depend_toymc: $(TOYMC:%.o=%.cxx) 
	rm -f ./.depend_toymc
	$(CXX) $(CPPFLAGS) $(ROOTCFLAGS) -MM $^ -MT $(TOYMC) > ./.depend_toymc

include .depend_toymc

# TOY MonteCarlo binary
###########################################################################

###########################################################################
# Generic development testing binary
tests/development.o: tests/development.cxx $(HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

bin/dev.app: tests/development.o
	@mkdir -p bin
	$(LD) $(LDFLAGS) ${ROOTGLIBS} ${ROOTLIBS} -o $@ $^

.PHONY: dev
dev:bin/dev.app
# Generic development testing binary
###########################################################################


############################################################################
# Generic partition testing binary
tests/partition.o: tests/partition.cxx $(HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

bin/part.app: tests/partition.o
	@mkdir -p bin
	$(LD) $(LDFLAGS) ${ROOTGLIBS} ${ROOTLIBS} -o $@ $^

.PHONY: part
part:bin/part.app
# Generic development testing binary
###########################################################################


.PHONY: clean
clean:
	@rm -f bin/*.app
	@rm -f tests/*.o
	@rm -f ToyMC/*.o
	@rm -f ToyMC/cint_dictionary.*
	@rm -f ToyMC/*.pcm
	@rm -f .depend_*


ToyMC/cint_dictionary.cxx: ToyMC/LinkDef.h
	rootcint -v4 -f $@ -c ToyMC/ToyMCEvent.h ToyMC/ToyMCParticle.h ToyMC/ToyMCGenerator.h $<
	-@[ -e "ToyMC/cint_dictionary_rdict.pcm" ] && mv -f ToyMC/cint_dictionary_rdict.pcm bin/		# ROOT 6

ToyMC/cint_dictionary.o: ToyMC/cint_dictionary.cxx
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(ROOTCFLAGS) $< -o $@

# include the toymc.o because it has the implementation of the logger
lib/ToyMC.so: ToyMC/cint_dictionary.o $(TOYMC)
	@mkdir -p lib
	$(LD) $(SOFLAGS) $(LDFLAGS) $(ROOTLIBS) $^ -o $@

.PHONY: rootlib
rootlib: lib/ToyMC.so
# DO NOT DELETE
