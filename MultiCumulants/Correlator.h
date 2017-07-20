#ifndef MULTICUMULANTS_CORRELATOR_H
#define MULTICUMULANTS_CORRELATOR_H

#include "MultiCumulants/QVectorSet.h"
#include "MultiCumulants/QTerms.h"

#include "vendor/loguru/loguru.hpp"

#include <set>
#include <utility> 

namespace cumulant{

    class Correlator{

    public:
        bool DEBUG = false;
        Complex v;
        Complex w;

        NativeMask _m;
        size_t _n;

        Correlator() : v(0, 0), w(0, 0) {

        }
        Correlator( NativeMask m, size_t n, QVectorMap &qvm) : v(0, 0), w(0, 0) {
            build( m, n, qvm );
        }

        int factorial( int n ){
            if ( n <= 1 ) return 1;
            return n * factorial( n - 1 );
        }

        void build( NativeMask m, size_t n, QVectorMap &qvm){
            // just save for printing
            _m = m;
            _n = n;
            LOG_F( INFO, "computing correlator for n=%zu", n );

            for ( auto qqq : qvm ){
                LOG_S(INFO) << qqq.first.to_string() << " : " << qqq.second.toString();
            }
            

            LOG_F( INFO, "m=%s", std::bitset<8>(m).to_string().c_str() );;

            auto lut = NativeMaskLUTs[ n-2 ];    
            auto coefLut = CoefficientKs[ n-2 ];
            size_t nTerms = lut.size();

            std::set< std::set< NativeMask > > visited;
        
            LOG_IF_F( INFO, DEBUG, "nTerms = %zu", nTerms );

            Complex qv(0, 0);
            Complex qw(0, 0);

            double totalK;
            // Loop over the number of terms in the correlator
            for ( size_t i = 0; i < nTerms; i++ ){
                LOG_F( INFO, "TERM %lu\n", i );
                Complex tv;
                Complex tw;
                std::string msg = "";
                std::string qmsg = "";
                // loop over the # of products in each term (maximum of n)
                size_t nVisited = 0;
                size_t nTest = 0;
                std::set<NativeMask> uTerms;
                for ( size_t j = 0; j < lut[ i ].size(); j++ ){
                    
                    NativeMask tm = maskAndCompactify( lut[i][j], m, 0, n );

                    if ( 0 == tm ) continue;
                    uTerms.insert( tm );
                    // nTest++;
                    // auto testPair = std::make_pair( lut[ i ].size(), tm );
                    std::bitset<MAX_SET_SIZE> bs( tm );
                    //     // LOG_F( INFO, " lut[ i ].size()=%lu, tm=%s", lut[ i ].size(), bs.to_string().c_str() );;
                    
                    // if (visited.count( testPair ) > 0 ) {
                    //     // LOG_F( INFO, "Visited"  );;
                        
                    //     // continue;
                    //     nVisited++;
                    // }
                    
                    
                    


                    

                    
                    // LOG_F( INFO, "bs = %s", bs.to_string().c_str() );
                    
                    if ( qvm.count( bs ) == 0 ){
                        LOG_F( INFO, "NOT FOUND" );
                    }

                    auto q = qvm[ bs ];
                    LOG_IF_F( INFO, DEBUG, "part tv=%f + i%f", q.getQV().real(), q.getQV().imag() );
                    LOG_IF_F( INFO, DEBUG, "part tw=%f + i%f", q.getW().real(), q.getW().imag() );

                    double ck = (pow(-1, q._i) * factorial(q._i));
                    // LOG_F( INFO, "nTerms=%lu, k=%zu, (-1)^k(k-1)! = %f * %s", nTerms, q._i+1, ck, bs.to_string().c_str() );


                    if ( 0 == j ){
                        totalK = ck;
                        msg = bs.to_string();
                        // qmsg = "(" + std::to_string( q.getQV().real() ) + ", " + std::to_string( q.getQV().imag() ) + ")";
                        qmsg = q.toString();
                        tv = q.getQV() ;
                        tw = q.getW() ;
                    }
                    else {
                        totalK *=ck;
                        msg += " " + bs.to_string();
                        // qmsg += "*(" + std::to_string( q.getQV().real() ) + ", " + std::to_string( q.getQV().imag() ) + ")";
                        qmsg += "* " + q.toString();
                        tv *= q.getQV() ;
                        tw *= q.getW() ;
                    }
                }
                
                
                if ( visited.count( uTerms ) > 0 ){
                    LOG_F( INFO, "should skip : %s", msg.c_str() );
                    continue;
                }
                visited.insert( uTerms );
                // LOG_IF_F( INFO, DEBUG, "coeff = %ld", k );
                // LOG_F( INFO, "visited? %lu == %s", visited.count( msg ), msg.c_str() );
                // visited.insert( msg );

                qv += tv * totalK;
                qw += tw  * totalK;

                // LOG_IF_F( INFO, DEBUG, "tv=%f + i%f", tv.real(), tv.imag() );
                // LOG_IF_F( INFO, DEBUG, "tw=%f + i%f", tw.real(), tw.imag() );

                LOG_F( INFO, "%f * %s", totalK, msg.c_str() );
                LOG_F( INFO, "%f * %s", totalK, qmsg.c_str() );
                


            }
            this->v = qv;
            this->w = qw;

            LOG_IF_F( INFO, DEBUG, "qv=%f + i%f", qv.real(), qv.imag() );
            LOG_IF_F( INFO, DEBUG, "qw=%f + i%f", qw.real(), qw.imag() );
            LOG_IF_F( INFO, DEBUG, "Finished computing n=%zu correlator", n );
        }

        Complex calculate(  ){
            return (this->v.real() / this->w.real());
        }


        inline NativeMask maskAndCompactify( NativeMask &im, NativeMask &mm, size_t start = 0, size_t stop = 8 ){

            // first do a quick check for validity
            
            NativeMask rm = im & mm;
            // require at least one bit to be 1
            if ( 0 == rm ) return 0;
            // then require that bits outside of mask are falsy
            if ( (im & (!mm)) > 0 ) return 0;

            NativeMask frm = 0;
            size_t n = 0;
            for ( size_t i = start; i < stop; i++ ){
                NativeMask ithbit = (1 << i);
                NativeMask nthbit = (1 << n);
                if ( ithbit & mm ){
                    // LOG_F( INFO, "Setting %luth bit", n );
                    if ( rm & ithbit )
                        frm |= (nthbit);
                    n++;
                }
            }

            // LOG_IF_F( INFO, DEBUG, "%s = (im=%s, mm=%s, rm=%s)", std::bitset<8>( frm ).to_string().c_str(), std::bitset<8>(im).to_string().c_str(), std::bitset<8>(mm).to_string().c_str(), std::bitset<8>(rm).to_string().c_str() );
            
            return frm;
        } // maskAndCompactify
    
        std::string toString(){
            std::string s = "";
            s += "Correlator( m=" + std::bitset<8>( _m ).to_string() + ", n=" + std::to_string(_n) + " )[ ";
            s += "v=" + std::to_string( v.real() ) + " + " + std::to_string( v.imag() ) + "i, ";
            s += "w=" + std::to_string( w.real() ) + " + " + std::to_string( w.imag() ) + "i";
            s += " ]";
            return s;

        }
    
    
    };
}


#endif
