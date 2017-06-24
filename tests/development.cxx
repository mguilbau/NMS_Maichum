#include "MultiCumulants/QVector.h"
#include "MultiCumulants/Result.h"
#include "MultiCumulants/Subsets.h"
#include "MultiCumulants/Algorithm.h"

#include "vendor/cmdline.h"

#include <iostream>
#include <cstdlib>
#include <string>
#include <random>
#include <chrono>
#include <fstream>

using namespace std;

void checkParam(int argc, char** argv);

void cumulants();
void mc();
void benchmark( int impl = 1, size_t horder = 4);

int 
main(int argc, char** argv) {

	cmdline::parser parser;

	parser.add("mc", '\0', "Run ToyMc");
	parser.add("cumulants", '\0', "Run Cumulants");
	parser.add<int>("benchmark", '\0', "Run Benchmarks", false, -1);
	parser.add<size_t>("horder", '\0', "harmonic order (benchmarks)", false, 4);

	parser.parse_check( argc, argv );


	loguru::set_thread_name("MAIN");
	// logs everything to the debug.log file every run
	loguru::add_file("bin/debug.log", loguru::Truncate, loguru::Verbosity_MAX);

	// sometimes the "stream" form of the logger is more convenient, use LOG_S( LEVEL ) << "MESSAGE";
	// No need for an endl at the end of a line
	
	if ( parser.exist( "mc" ) ){
		mc();
	} else if ( parser.exist( "cumulants" ) ){
		cumulants();
	} else if ( parser.get<int>( "benchmark" ) >= 1 ){
		// benchmark( parser.get<int>( "benchmark" ), parser.get<size_t>( "horder" ) );
	} else {
		cout << parser.usage() << endl;
	}
	return 0;
}

void cumulants(){
	LOG_F( INFO, "Cumulants" );

	const size_t order = 3;

	cumulant::Subset s1(2);
	s1.set(0, "pT", 0., 10.0);
	s1.set(1, "eta",0,10);
	cumulant::Subset s2(2);
	s2.set(0, "pT", 0., 10.0);
	s2.set(1, "eta",-10, 0);
	cumulant::Subset s3(2);
	s3.set(0, "pT", 0., 10.0);
	s3.set(1, "eta",-10, 10);
	cumulant::Subset s4(2);
	s4.set(0, "pT", 0., 10.0);
	s4.set(1, "eta",-10, 10);
	LOG_S(INFO) << s1.toString();
	LOG_S(INFO) << s2.toString();
	LOG_S(INFO) << s3.toString();
	LOG_S(INFO) << s4.toString();

	cumulant::Set qvset(order);
	qvset.setSubsetParams(0,s1);
	qvset.setSubsetParams(1,s2);
	qvset.setSubsetParams(2,s3);
	qvset.setSubsetParams(3,s4);
	LOG_S(INFO) << qvset.toString();


	double PI = 3.14159265358;
	std::vector<double> phi;
	std::vector<double> eta;
	std::vector<double> pt;
	std::vector<double> w;

	std::random_device rd;
	std::mt19937 ephi(rd());
	std::mt19937 eeta(rd());
	std::mt19937 ept(rd());
	std::mt19937 ew(rd());
	std::uniform_real_distribution<> distphi(0, 2*PI);
	std::uniform_real_distribution<> disteta(-10, 10);
	std::uniform_real_distribution<> distpt(0, 10);
	std::uniform_real_distribution<> distw(0, 1);

	size_t NEvt = 100000;
	for (int n = 0; n < NEvt; ++n) {
		phi.push_back(distphi(ephi));
		eta.push_back(disteta(eeta));
		pt.push_back(distpt(ept));
		w.push_back(1.);
	}

	std::vector< vector<double> > val(NEvt, std::vector<double>(2,0.));
	for (int n = 0; n < NEvt; ++n) {
		val[n][0] = pt[n];
		val[n][1] = eta[n];
	}


	HarmonicVector h(order);
	h[0] = 2;
	h[1] = 2;
	h[2] = -2;
	h[3] = -2;

	cumulant::QVectorSet qv(h,qvset,false);
	qv.reset();

	LOG_S(INFO) << "\n" << qv.maskString() << endl;
	for (int n = 0; n < NEvt; ++n) {
	   // qv.generateMask(val[n]);
	   qv.fill(val[n], phi[n], w[n]);
	}
	// 
	
	LOG_S(INFO) << qv.print();
	qv.reset();

	cumulant::QVectorVector q = qv.getQ();

	cumulant::QTerms qt;

	qt.test( 4 );

}

void mc(){
	LOG_F( INFO, "ToyMc" );
}


//
void checkParam(int argc, char** argv)
{
		LOG_S(INFO) << "Number of parameters: " << argc;
		for(int ip=0; ip<argc; ++ip) 
		  LOG_S(INFO) << "Argument " << ip << ": " << argv[ip];
}