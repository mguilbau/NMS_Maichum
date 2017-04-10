// logging library
#define LOGURU_IMPLEMENTATION 1
#include "vendor/loguru/loguru.hpp"

#include <iostream>
#include <cstdlib>
#include <string>
#include <random>

#include "ToyMC/ToyMCGenerator.h"
#include "ToyMC/BranchReader.h"
#include "ToyMC/TClonesArrayReader.h"

#include <TFile.h>
#include <TDirectory.h>
#include <TH1.h>
#include <TSpline.h>
#include <TF1.h>

using namespace std;

void checkParam(int argc, char** argv);
void createDistributions();

int 
main(int argc, char** argv) {

	checkParam(argc, argv);

	loguru::set_thread_name("MAIN");
	// logs everything to the debug.log file every run
	loguru::add_file("bin/debug.log", loguru::Truncate, loguru::Verbosity_MAX);

	// sometimes the "stream" form of the logger is more convenient, use LOG_S( LEVEL ) << "MESSAGE";
	// No need for an endl at the end of a line
	
	//toymc::ToyMCGenerator g_const;
	//LOG_S(INFO) << g_const.toString();
	//g_const.generate(10000, "../output/ToyMCTTree_noVn_const.root");

	//toymc::ToyMCGenerator g_exp("PbPb", "exp", "");
	//LOG_S(INFO) << g_exp.toString();
	//g_exp.generate(10000, "../output/ToyMCTTree_noVn_exp.root");

	toymc::ToyMCGenerator g_const2("PbPb", "exp", "exp");
        g_const2.setRanges(0.3,3.0,-2.4,2.4,1,500);
	LOG_S(INFO) << g_const2.toString();
	g_const2.generate(10000000, "/Volumes/Elements/ToyMCdata/datafiles/cumulant/v1/ToyMCTTree_expVn_exp.root");

	//toymc::ToyMCGenerator g_const3("PbPb", "const", "exp");
	//LOG_S(INFO) << g_const3.toString();
	//g_const3.generate(10000, "../output/ToyMCTTree_expVn_const.root");


//        TFile* fin = TFile::Open("/Volumes/Elements/ToyMCdata/datafiles/cumulant/v1/ToyMCTTree_noVn_const.root");
//        TTree* tr = (TTree*) fin->Get("ToyMCTree");
//        BranchReader<toymc::ToyMCEvent> eventReader;
//        TClonesArrayReader<toymc::ToyMCParticle> plcsReader;
//        eventReader.setup(tr, "Event");
//        plcsReader.setup(tr, "Particles");
//
//        toymc::ToyMCEvent event;
//        toymc::ToyMCParticle plc;
//
//        int nentries = tr->GetEntries();
//        tr->Print();
//        int ievt = 0;
//        while(ievt < nentries)
//        {
//           plc.reset();
//
//           tr->GetEntry(ievt);
//           std::cout <<
//           "\rtoymc.cxx::INFO:: ievt = " << ievt
//           <<
//           " ~~~> " << std::setprecision(3) << (double)((double)ievt / (double)nentries)*100  << " %"
//           << std::flush;
//
//           int npart = plcsReader.N();
//           for(int ip = 0; ip < npart; ++ip)
//           {
//              plc = *plcsReader.get(ip);
//              //std::cout << plc.getpt() << std::endl;
//           }
//
//           ++ievt;
//        }
//        std::cout << std::endl;


	return 0;
}
//
void checkParam(int argc, char** argv)
{
        LOG_S(INFO) << "Number of parameters: " << argc;
        for(int ip=0; ip<argc; ++ip) 
          LOG_S(INFO) << "Argument " << ip << ": " << argv[ip];
}
