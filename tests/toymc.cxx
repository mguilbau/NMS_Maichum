#include "ToyMC/ToyMCGenerator.h"

#include <iostream>
#include <cstdlib>
#include <string>
#include <random>

// logging library
#define LOGURU_IMPLEMENTATION 1
#include "vendor/loguru/loguru.hpp"

using namespace std;

void checkParam(int argc, char** argv);

int 
main(int argc, char** argv) {

	checkParam(argc, argv);

	loguru::set_thread_name("MAIN");
	// logs everything to the debug.log file every run
	loguru::add_file("bin/debug.log", loguru::Truncate, loguru::Verbosity_MAX);

	// sometimes the "stream" form of the logger is more convenient, use LOG_S( LEVEL ) << "MESSAGE";
	// No need for an endl at the end of a line
	
	toymc::ToyMCGenerator g;
	LOG_S(INFO) << g.toString();
	g.generate(1000);

	return 0;
}
//
void checkParam(int argc, char** argv)
{
        LOG_S(INFO) << "Number of parameters: " << argc;
        for(int ip=0; ip<argc; ++ip) 
          LOG_S(INFO) << "Argument " << ip << ": " << argv[ip];
}
