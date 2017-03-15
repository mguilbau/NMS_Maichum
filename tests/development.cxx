#include "MultiCumulants/QVector.h"
#include "MultiCumulants/Subsets.h"

#include <iostream>
using namespace std;




int 
main(int argc, char** argv) {

	//loguru::set_thread_name("MAIN");
	// logs everything to the debug.log file every run
	//loguru::add_file("bin/debug.log", loguru::Truncate, loguru::Verbosity_MAX);


	QVector qv;
	Subset er;


	// sometimes the "stream" form of the logger is more convenient, use LOG_S( LEVEL ) << "MESSAGE";
	// No need for an endl at the end of a line
	//LOG_S(INFO) << qv.toString();
	//LOG_S(INFO) << er.toString();
	


	return 0;
}
