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

        Correlator() : v(0, 0), w(0, 0) {

        }
        Correlator( NativeMask m, QVectorMap &qvm) : v(0, 0), w(0, 0) {
            build( m, qvm );
        }

        int factorial( int n ){
            if ( n <= 1 ) return 1;
            return n * factorial( n - 1 );
        }

        size_t countSetBits(NativeMask m){
            size_t count = 0;
            while (m){
                m &= (m-1) ;
                count++;
            }
            return count;
        }

        // im is the one in the compressed space
        // mm is the mask in the full space
        // example :
        // given im=0101, and mm=00110101
        // returns rm = 00010001, ie 1st and 3rd set bits in mask
        NativeMask expandMask( NativeMask im, NativeMask mm, size_t start = 0, size_t stop = 8 ){

            vector<size_t> mlut;
            size_t n = 0;
            NativeMask rm = 0;
            for ( size_t i = start; i <stop; i++ ){
                NativeMask ithbit = (1 << i);
                NativeMask nthbit = (1 << n);
                if ( ithbit & mm ){
                    if ( nthbit & im )
                        rm |= ithbit;
                    n++; 
                }
            }
            return rm;
        }


        void build( NativeMask m, QVectorMap &qvm){
            // just save for printing
            _m = m;

            auto bm = std::bitset<8>( m );
            size_t maskBitsSet = countSetBits( m );
            LOG_IF_F( INFO, DEBUG, "nSetBits(mask) = %lu", maskBitsSet );

            auto lut = NativeMaskLUTs[ maskBitsSet-2 ];    

            size_t nTerms = lut.size();
            Complex qv(0, 0);
            Complex qw(0, 0);
            for ( size_t i = 0; i < nTerms; i++ ){
                LOG_IF_F( INFO, DEBUG, "\n\nTERM %lu", i );
                
                double totalK = 1.0;
                Complex tv(0,0);
                Complex tw(0,0);
                std::string msg = "";
                for ( size_t j = 0; j < lut[ i ].size(); j++ ){
                    NativeMask tm = lut[ i ][ j ];
                    NativeMask em = expandMask( tm, m );
                    
                    auto btm = std::bitset<8>( tm );
                    auto bem = std::bitset<8>( em );
                    
                    LOG_IF_F( INFO, DEBUG, "NativeMask=%s", std::bitset<8>( tm ).to_string().c_str() );
                    LOG_IF_F( INFO, DEBUG, "expandMask( im=%s, mm=%s ) = %s", btm.to_string().c_str(), bm.to_string().c_str(), bem.to_string().c_str() );

                    auto q = qvm[ bem ];
                    double ck = (pow(-1, q._i) * factorial(q._i));
                    totalK *= ck;
                    msg += " " + bem.to_string();
                    
                    if ( 0 == j ){
                        tv = q.getQV();
                        tw = q.getW();
                    } else {
                        tv *= q.getQV();
                        tw *= q.getW();
                    }
                    

                } // loop on j

                LOG_IF_F( INFO, DEBUG, "%f * %s", totalK, msg.c_str() );

                qv += tv * totalK;
                qw += tw * totalK;

            } // loop on i terms
    
            this->v = qv;
            this->w = qw;

            return;

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
            s += "Correlator( m=" + std::bitset<8>( _m ).to_string() + " )[ ";
            s += "v=" + std::to_string( v.real() ) + " + " + std::to_string( v.imag() ) + "i, ";
            s += "w=" + std::to_string( w.real() ) + " + " + std::to_string( w.imag() ) + "i";
            s += " ]";
            return s;

        }
    
    
    };
}


#endif
