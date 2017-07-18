#include "MultiCumulants/QVector.h"
#include "MultiCumulants/QVectorSet.h"
#include "MultiCumulants/QTerms.h"
#include "MultiCumulants/Result.h"
#include "MultiCumulants/Subsets.h"
#include "MultiCumulants/Algorithm.h"
#include "MultiCumulants/Correlator.h"

#include "vendor/cmdline.h"

#include <iostream>
#include <cstdlib>
#include <string>
#include <random>
#include <chrono>
#include <fstream>
#include <iterator>
 
using namespace std;

void checkParam(int argc, char** argv);
void cumulants( size_t h );

#define LAST(k,n) ((k) & ((1<<(n))-1))
#define MID(k,m,n) LAST((k)>>(m),((n)-(m)))

NativeMask maskAndCompactify( NativeMask im, NativeMask mm ){

	// LOG_F( INFO, "(im=%s, mm=%s)", std::bitset<8>(im).to_string().c_str(), std::bitset<8>(mm).to_string().c_str() );
	NativeMask rm = im & mm;
	// LOG_F( INFO, "(im=%s, mm=%s, rm=%s)", std::bitset<8>(im).to_string().c_str(), std::bitset<8>(mm).to_string().c_str(), std::bitset<8>(rm).to_string().c_str() );
	NativeMask frm = 0;
	size_t n = 0;
	for ( size_t i = 0; i < 8; i++ ){
		NativeMask ithbit = (1 << i);
		NativeMask nthbit = (1 << n);
		// LOG_F( INFO, "Testing ith bit : %s", std::bitset<8>( ithbit ).to_string().c_str() );
		// LOG_F( INFO, "nth bit : %s", std::bitset<8>( nthbit ).to_string().c_str() );
		if ( ithbit & mm ){
			// LOG_F( INFO, "Setting %luth bit", n );
			if ( rm & ithbit )
				frm |= (nthbit);
			n++;
		}
	}

	LOG_F( INFO, "%s = (im=%s, mm=%s, rm=%s)", std::bitset<8>( frm ).to_string().c_str(), std::bitset<8>(im).to_string().c_str(), std::bitset<8>(mm).to_string().c_str(), std::bitset<8>(rm).to_string().c_str() );
	// LOG_F( INFO, "result = %s", std::bitset<8>( frm ).to_string().c_str() );
	return frm;
}


int 
main(int argc, char** argv) {

	cmdline::parser parser;
	parser.add<size_t>("h", '\0', "harmonic order (benchmarks)", false, 4);
	parser.parse_check( argc, argv );
	
	loguru::set_thread_name("MAIN");
	


	// logs everything to the debug.log file every run
	loguru::add_file("bin/debug.log", loguru::Truncate, loguru::Verbosity_MAX);
	
	cumulants( parser.get<size_t>( "h" ) );
	
	return 0;
}

void cumulants( size_t order ){
	LOG_F( INFO, "Cumulants( order=%zu )", order  );

	cumulant::Set qvset(order);
	for ( size_t i = 0; i < order; i++ ){
		cumulant::Subset ss(2);
		ss.set(0, "pT", 0., 10.0);
		ss.set(1, "eta",0,10);

		qvset.setSubsetParams( i, ss );
		// LOG_S(INFO) << ss.toString();
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

	size_t NEvt = 100000;
	for ( size_t n = 0; n < NEvt; ++n) {
		phi.push_back(distphi(ephi));
		eta.push_back(disteta(eeta));
		pt.push_back(distpt(ept));
		w.push_back(1.);
	}

	std::vector< vector<double> > val(NEvt, std::vector<double>(2,0.));
	for ( size_t n = 0; n < NEvt; ++n) {
		val[n][0] = pt[n];
		val[n][1] = eta[n];
	}


	LOG_F( INFO, "Creating HarmonicVector(%zu)", order );
	HarmonicVector h(order);
	h[0] = 2;
	h[1] = 2;
	h[2] = -2;
	h[3] = -2;
	LOG_F( INFO, "HarmonicVector.size()=%zu", h.size() );


	cumulant::QVectorSet qv(h,qvset,false);
	qv.reset();

	LOG_S(INFO) << "\n" << qv.maskString() << endl;
	for (size_t n = 0; n < NEvt; ++n) {
	   // qv.generateMask(val[n]);
	   qv.fill(val[n], phi[n], w[n]);
	}
	// 
	
	LOG_S(INFO) << qv.print();
	// qv.reset();

	cumulant::QVectorMap& q = qv.getQ();

	// cumulant::QTerms qt;
	// qt.generate( order, true );

	cumulant::Correlator C( 0b111010, order, q );


	// size_t nTerms = sizeof( cumulant::QTERMS_h8 ) / sizeof( cumulant::QTERMS_h8[0] );
	// LOG_F( INFO, "computing 8th order correlator" );
	// long int tc = 1;
	// for ( size_t i = 0; i < nTerms; i++ ){
	// 	tc = cumulant::QTERMS_h8_KCOEFF[i];
	// 	NativeMask &nm = cumulant::QTERMS_h8[i][0];
	// 	std::bitset<MAX_SET_SIZE> bs( nm );
	// 	NativeMask nm2 = bs.to_ullong();
	// 	if ( nm != nm2 )
	// 		LOG_F( INFO, "mismatch" );

	// 	// LOG_F( INFO, "kcoeff=%ld", cumulant::QTERMS_h8_KCOEFF[i] );
	// }
	// LOG_F( INFO, "finished computing 8th order correlator, tc=%ld", tc );

	// LOG_F( INFO, "sizeof(TERMSMAP 8) = %lu", sizeof( cumulant::QTERMS_h8 ) / sizeof( cumulant::QTERMS_h8[0] ) );
	// LOG_F( INFO, "%s", std::bitset<MAX_SET_SIZE>( cumulant::QTERMS_h4[0][0]).to_string().c_str() );

	// int a = 0xabcdef;
    // printf("%x\n",  MID(a,4,16));
	// LOG_F( INFO, "%s ==> %s", std::bitset<16>( a ).to_string().c_str(), std::bitset<16>( MID(a,4,16) ).to_string().c_str() );


	// NativeMask rm = maskAndCompactify( 0b10101010, 0b00001010 );
	// rm = maskAndCompactify( 0b10101010, 0b00001011 );
	// rm = maskAndCompactify( 0b10101010, 0b00001110 );
	// rm = maskAndCompactify( 0b10101010, 0b00111010 );
}

void checkParam(int argc, char** argv)
{
		LOG_S(INFO) << "Number of parameters: " << argc;
		for(int ip=0; ip<argc; ++ip) 
		  LOG_S(INFO) << "Argument " << ip << ": " << argv[ip];
}
