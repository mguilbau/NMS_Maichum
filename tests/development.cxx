#include "MultiCumulants/QVector.h"
#include "MultiCumulants/QVectorSet.h"
#include "MultiCumulants/QTerms.h"
#include "MultiCumulants/Result.h"
#include "MultiCumulants/Subsets.h"
#include "MultiCumulants/Algorithm.h"
#include "MultiCumulants/Correlator.h"
#include "MultiCumulants/Cumulant.h"

#include <correlations/Types.hh>
#include <correlations/Result.hh>
#include <correlations/QVector.hh>
#include <correlations/recursive/FromQVector.hh>
#include <correlations/recurrence/FromQVector.hh>

#include "vendor/cmdline.h"

#include <iostream>
#include <cstdlib>
#include <string>
#include <random>
#include <chrono>
#include <fstream>
#include <iterator>
 
using namespace std;

class Events {
  public:
    Events( )             { setMult( 500 ); }
    Events( size_t mult ) { setMult( mult ); }
    ~Events() {};

    void setMult( size_t mult ) 
    {
       mult_ = mult;
       phi_.resize(mult);
       pt_.resize(mult);
       eta_.resize(mult);
       w_.resize(mult);
    } 
    size_t getMult() { return mult_; }

    void setPhi( size_t m, double phi ) 
    {  
       phi_[m] = phi;
    }
    void setPt( size_t m, double pt )  
    {  

       pt_[m] = pt;
    }
    void setEta( size_t m, double eta ) 
    {  
   
       eta_[m] = eta;
    } 
    void setW( size_t m, double w )   
    {  
       w_[m] = w;
    }
    std::vector<double> getPhi() { return phi_; }
    std::vector<double> getPt()  { return pt_; }
    std::vector<double> getEta() { return eta_; } 
    std::vector<double> getW()   { return w_; }

  private:
    size_t mult_;
    std::vector<double> phi_;
    std::vector<double> eta_;
    std::vector<double> pt_;
    std::vector<double> w_;

};

void checkParam(int argc, char** argv);
void randomfunction( size_t h );
void cumulants( size_t order, size_t nset, size_t nevents, size_t mult, int harm );

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
	parser.add<size_t>("order", '\0', "cumulant order (benchmarks)", false, 4);
	parser.add<size_t>("nset", '\0', "number of subsets", false, 4);
	parser.add<size_t>("nevents", '\0', "number of events", false, 1);
	parser.add<size_t>("multmax", '\0', "event multiplicity", false, 1);
	parser.add<int>("harm", '\0', "harmonic order", false, 2);
	parser.parse_check( argc, argv );
	
	loguru::set_thread_name("MAIN");
	

	// logs everything to the debug.log file every run
	loguru::add_file("bin/debug.log", loguru::Truncate, loguru::Verbosity_MAX);
	
	cumulants( parser.get<size_t>( "order" ),
                   parser.get<size_t>( "nset" ),
                   parser.get<size_t>( "nevents" ),
                   parser.get<size_t>( "multmax" ),
                   parser.get<int>( "harm" ) );
	
	return 0;
}

void cumulants( const size_t order, size_t nset, size_t nevents,  size_t mult, int harm ) 
{
	LOG_F( INFO, "Cumulants( order=%zu )", order  );
	cumulant::Set set(nset);
	for ( size_t i = 0; i < nset; i++ ) {
		cumulant::Subset ss(2);
		ss.set(0, "pT",  0.3, 3.0);
		ss.set(1, "eta",-2.4, 2.4);

		set.setSubsetParams( i, ss );
		LOG_S(INFO) << ss.toString();
	}
	LOG_S(INFO) << set.toString();

        //Get random distributions
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
	std::uniform_real_distribution<> distpt(0.3, 3.);
	std::uniform_real_distribution<> disteta(-2.4, 2.4);
	std::uniform_real_distribution<> distw(0, 1);

        //Fill Events
	size_t nEvt  = nevents;
	size_t nMult = mult;

        std::vector<Events> evtvec( nEvt, Events( nMult ) );

	LOG_S(INFO) << "Looping over " << nEvt << " with mult = " << nMult;
	for ( size_t n = 0; n < nEvt; ++n) {
	   for ( size_t m = 0; m < nMult; ++m) {
		evtvec[n].setPhi( m, distphi(ephi) );
		evtvec[n].setPt(  m, distpt(ept) );
		evtvec[n].setEta( m, disteta(eeta) );
		evtvec[n].setW(   m, distw(ew) );
		//evtvec[n].setW(   m, 0. );
		//evtvec[n].setW(   m, 0.1 );
		//evtvec[n].setW(   m, 1. );
           }
        }

	LOG_F( INFO, "Creating HarmonicVector(%zu)", order );
	HarmonicVector h(order);
	for ( size_t hh = 0; hh < order; ++hh) {
           if( hh < order/2 ) h[hh] =  1*harm;
           else               h[hh] = -1*harm;
        }
	LOG_F( INFO, "HarmonicVector.size()=%zu", h.size() );

	std::vector<double> val(2,0.);
        std::bitset<8> bitset;

        //Init standard method
        correlations::QVector qstd(0, 0, true);
        correlations::HarmonicVector hstd;
        correlations::FromQVector *cstd;
        hstd = correlations::HarmonicVector(8);
        hstd[0] =  1*harm;
        hstd[1] = -1*harm;
        hstd[2] =  1*harm;
        hstd[3] = -1*harm;
        hstd[4] =  1*harm;
        hstd[5] = -1*harm;
        hstd[6] =  1*harm;
        hstd[7] = -1*harm;
        qstd.resize(hstd);
        cstd = new correlations::recurrence::FromQVector(qstd);
        correlations::Result rN2;
        correlations::Result rN4;
        correlations::Result rN6;
        correlations::Result rN8;

        //cumulant building
	cumulant::QVectorSet qv(h, set, true);

	for (size_t n = 0; n < nEvt; ++n) {
	   //LOG_S(INFO) << "\n" << qv.maskString() << endl;
	   qstd.reset();
	   qv.reset();

           LOG_S(INFO) << "Event number = " << n+1;
	   for (size_t m = 0; m < nMult; ++m) {
	       val[0] = evtvec[n].getPt()[m];
	       val[1] = evtvec[n].getEta()[m];
               
               LOG_S(INFO) << "~~~~> Track number = " << m+1 << ", phi = " << evtvec[n].getPhi()[m] 
                                                             << ", pt  = " << evtvec[n].getPt()[m]
                                                             << ", eta = " << evtvec[n].getEta()[m]
                                                             << ", w   = " << evtvec[n].getW()[m];

               //Cumulate
               qstd.fill( evtvec[n].getPhi()[m], evtvec[n].getW()[m] );
	       qv.fill( val, evtvec[n].getPhi()[m], evtvec[n].getW()[m] );
	   }
	   LOG_S(INFO) << qv.print();
	   cumulant::QVectorMap& q = qv.getQ();
	   std::vector<cumulant::Correlator> cn;
           cn.resize(order/2);

           //c2
           bitset.set(0);
           bitset.set(4);
           LOG_S(INFO) << " bitset cn2 = " << bitset;
           cn[0] = cumulant::Correlator( bitset.to_ulong(), q );
           //Bilandzic code
           rN2 = cstd->calculate(2, hstd);
           LOG_S(INFO) << " New FWK ~~> c_{" << harm << "}{2} = " << cn[0].calculate() << endl;
           LOG_S(INFO) << " Old FWK ~~> c_{" << harm << "}{2} = " << rN2.eval() << endl;
           //c4
           if(order > 2)
           {
              bitset.set(1);
              bitset.set(5);
              LOG_S(INFO) << " bitset cn4 = " << bitset;
              cn[1] = cumulant::Correlator( bitset.to_ulong(), q );
              //Bilandzic code
              rN4 = cstd->calculate(4, hstd);
              LOG_S(INFO) << " New FWK ~~> c_{" << harm << "}{4} = " << cn[1].calculate() << endl;
              LOG_S(INFO) << " Old FWK ~~> c_{" << harm << "}{4} = " << rN4.eval() << endl;
           }
           //c6
           if(order > 4)
           {
              bitset.set(2);
              bitset.set(6);
              LOG_S(INFO) << " bitset cn6 = " << bitset;
              cn[2] = cumulant::Correlator( bitset.to_ulong(), q );
              //Bilandzic code
              rN6 = cstd->calculate(6, hstd);
              LOG_S(INFO) << " New FWK ~~> c_{" << harm << "}{6} = " << cn[2].calculate() << endl;
              LOG_S(INFO) << " Old FWK ~~> c_{" << harm << "}{6} = " << rN6.eval() << endl;
           }
           //c8
           if(order > 6)
           {
              bitset.set(3);
              bitset.set(7);
              LOG_S(INFO) << " bitset cn8 = " << bitset;
              cn[3] = cumulant::Correlator( bitset.to_ulong(), q );
              //Bilandzic code
              rN8 = cstd->calculate(8, hstd);
              LOG_S(INFO) << " New FWK ~~> c_{" << harm << "}{8} = " << cn[3].calculate() << endl;
              LOG_S(INFO) << " Old FWK ~~> c_{" << harm << "}{8} = " << rN8.eval() << endl;
           }

           bitset.reset();
           LOG_S(INFO) << " bitset after reset = " << bitset << endl;
	}
}


void randomefunction( size_t order ){
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

	cumulant::Correlator C( 58, q );

	cumulant::Cumulant cum( 7 );
	cum.buildCorrelators( q );

	cum.buildCumulant();


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
