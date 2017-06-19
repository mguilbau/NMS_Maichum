#include "MultiCumulants/QVector.h"
#include "MultiCumulants/Result.h"
#include "MultiCumulants/Subsets.h"
#include "MultiCumulants/Algorithm.h"
#include "MultiCumulants/Correlator.h"

#include <correlations/Types.hh>
#include <correlations/Result.hh>
#include <correlations/QVector.hh>
#include <correlations/recursive/FromQVector.hh>
#include <correlations/recurrence/FromQVector.hh>

#include "ToyMC/ToyMCGenerator.h"
#include "ToyMC/BranchReader.h"
#include "ToyMC/TClonesArrayReader.h"
#include "ToyMC/ToyMCEvent.h"
#include "ToyMC/ToyMCParticle.h"
#include "ToyMC/BranchWriter.h"
#include "ToyMC/TClonesArrayWriter.h"

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
void mc(int harm,
        std::string system,  std::string partDist, std::string vnDist,
        double ptMin,        double ptMax,
        double etaMin,       double etaMax,
        double multMin,      double multMax,
        int    nEvt,         bool   vnFluct,
        std::string version, std::string outFileName, int nevts);
void benchmark( int impl = 1, size_t horder = 4);

int 
main(int argc, char** argv) {

	cmdline::parser parser;

	parser.add("mc", '\0', "Run ToyMc");
	parser.add("cumulants", '\0', "Run Cumulants");
	parser.add<int>("benchmark", '\0', "Run Benchmarks", false, -1);
	parser.add<size_t>("horder", '\0', "harmonic order (benchmarks)", false, 4);
	parser.add<int>("nevent", '\0', "number of generated mc events", false, 1);

	parser.parse_check( argc, argv );


	loguru::set_thread_name("MAIN");
	// logs everything to the debug.log file every run
	loguru::add_file("bin/debug.log", loguru::Truncate, loguru::Verbosity_MAX);

	// sometimes the "stream" form of the logger is more convenient, use LOG_S( LEVEL ) << "MESSAGE";
	// No need for an endl at the end of a line
	
	if ( parser.exist( "mc" ) ){
		mc( 2,
                    "PbPb", "const", "const", 
                    0.3, 3.0, 
                    -2.4, 2.4, 
                    10, 40, 
                    10000, true, "v6", "10k_test_genAnalyze", parser.get<int>( "nevent" ) );
	} else if ( parser.exist( "cumulants" ) ){
		cumulants();
	} else if ( parser.get<int>( "benchmark" ) >= 1 ){
		benchmark( parser.get<int>( "benchmark" ), parser.get<size_t>( "horder" ) );
	} else {
		cout << parser.usage() << endl;
	}
	return 0;
}

void 
cumulants(){
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

	cumulant::impl2::QVectorSet qv(h,qvset,false);
	qv.reset();

	LOG_S(INFO) << "\n" << qv.maskString() << endl;
	for (int n = 0; n < NEvt; ++n) {
	   // qv.generateMask(val[n]);
	   qv.fill(val[n], phi[n], w[n]);
	}
	// 
	
	LOG_S(INFO) << qv.print();
	qv.reset();

	cumulant::impl2::QVectorVector q = qv.getQ();

}

void mc(int harm,
        std::string system,  std::string partDist, std::string vnDist,
        double ptMin,        double ptMax,
        double etaMin,       double etaMax,
        double multMin,      double multMax,
        int    nEvt,         bool   vnFluct,
        std::string version, std::string outFileName, int nevts = -1)
{
	toymc::ToyMCGenerator g(system, partDist, vnDist);
        g.setRanges(ptMin, ptMax, etaMin, etaMax, multMin, multMax);
	LOG_S(INFO) << g.toString();
	if(vnFluct) g.setFlowFluctuations();
}


//
void checkParam(int argc, char** argv)
{
	LOG_S(INFO) << "Number of parameters: " << argc;
	for(int ip=0; ip<argc; ++ip) 
	LOG_S(INFO) << "Argument " << ip << ": " << argv[ip];
}



void benchmark( int impl, size_t order ){
	LOG_F( INFO, "BENCHMARK on filling QVectors" );


	cumulant::Set qvset(order);

	for ( size_t i = 0; i < order; i++ ){
		cumulant::Subset s(2);
		s.set(0, "pT", 0., 10.0);
		s.set(1, "eta",0,10);
		LOG_S(INFO) << s.toString();
		qvset.setSubsetParams(i,s);
	}
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


	size_t NEvents = 100;
	size_t NTrk = 10000;

	for (size_t n = 0; n < NTrk; ++n) {
		phi.push_back(distphi(ephi));
		eta.push_back(disteta(eeta));
		pt.push_back(distpt(ept));
		w.push_back(1.);
	}

	std::vector< vector<double> > val(NTrk, std::vector<double>(2,0.));
	for (size_t n = 0; n < NTrk; ++n) {
		val[n][0] = pt[n];
		val[n][1] = eta[n];
	}


	HarmonicVector h(order);
	h[0] = 2;
	h[1] = 2;
	h[2] = -2;
	h[3] = -2;

	

	stringstream name;
	name << "bench_" << impl << "_" << order << ".dat";
	ofstream benchFile( name.str().c_str() );
	benchFile << "implementation=" << impl << "\norder=" << order << "\n";
	benchFile << "NTrk/Event=" << NTrk << "\n";

	using namespace std::chrono;
	
	LOG_F( INFO, "%d\n%lu", impl, order );

	if ( 1 == impl ){
		cumulant::impl1::QVectorSet qv1(h,qvset,false);
		qv1.reset();

		for ( size_t iEvent = 0 ; iEvent < NEvents; iEvent++ ){
			high_resolution_clock::time_point t1 = high_resolution_clock::now();
			for (size_t n = 0; n < NTrk; ++n) {
				qv1.fill(val[n], phi[n], w[n]);
			}
			high_resolution_clock::time_point t2 = high_resolution_clock::now();
			duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
			benchFile << "duration=" << time_span.count() << "\n";
			LOG_F( INFO, "duration(%lu Tracks)=%f", NTrk, time_span.count() );
		}
	}

	if ( 2 == impl ){
		cumulant::impl2::QVectorSet qv2(h,qvset,false);
		qv2.reset();

		LOG_S( INFO ) << "\n" << qv2.maskString() << endl;

		for ( size_t iEvent = 0 ; iEvent < NEvents; iEvent++ ){
			high_resolution_clock::time_point t1 = high_resolution_clock::now();
			for (size_t n = 0; n < NTrk; ++n) {
				qv2.fill(val[n], phi[n], w[n]);
			}
			high_resolution_clock::time_point t2 = high_resolution_clock::now();
			duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
			LOG_F( INFO, "duration=%f", time_span.count() );
			benchFile << "duration=" << time_span.count() << "\n";
		}
	}
	

	if ( 3 == impl ){
		cumulant::impl3::QVectorSet qv3(h,qvset,false);
		qv3.reset();

		for ( size_t iEvent = 0 ; iEvent < NEvents; iEvent++ ){
			high_resolution_clock::time_point t1 = high_resolution_clock::now();
			for (size_t n = 0; n < NTrk; ++n) {
				qv3.fill(val[n], phi[n], w[n]);
			}
			high_resolution_clock::time_point t2 = high_resolution_clock::now();
			duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
			LOG_F( INFO, "duration=%f", time_span.count() );
			benchFile << "duration=" << time_span.count() << "\n";
		}
	}

	if ( 4 == impl ){
		cumulant::impl4::QVectorSet qv4(h,qvset,false);
		qv4.reset();

		for ( size_t iEvent = 0 ; iEvent < NEvents; iEvent++ ){
			high_resolution_clock::time_point t1 = high_resolution_clock::now();
			for (size_t n = 0; n < NTrk; ++n) {
				qv4.fill(val[n], phi[n], w[n]);
			}
			high_resolution_clock::time_point t2 = high_resolution_clock::now();
			duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
			LOG_F( INFO, "duration=%f", time_span.count() );
			benchFile << "duration=" << time_span.count() << "\n";
		}
	}
	

	benchFile.close();
	// LOG_S(INFO) << "\n" << qv.maskString() << endl;


	

	

}
