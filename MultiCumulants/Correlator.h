#ifndef MULTICUMULANTS_CORRELATOR_H
#define MULTICUMULANTS_CORRELATOR_H

#include "MultiCumulants/QVectorSet.h"
#include "MultiCumulants/QTerms.h"

#include "vendor/loguru/loguru.hpp"

namespace cumulant{

    class Correlator{

    public:
        bool DEBUG = false;
        Complex v;
        Complex w;
        Correlator() : v(0, 0) {

        }
        Correlator( NativeMask m, size_t n, QVectorMap &qvm) : v(0, 0) {
            build( m, n, qvm );
        }

        void build( NativeMask m, size_t n, QVectorMap &qvm){
            LOG_F( INFO, "computing correlator for n=%zu", n );

            auto lut = NativeMaskLUTs[ n-2 ];    
            size_t nTerms = lut.size();
        
            LOG_IF_F( INFO, DEBUG, "nTerms = %zu", nTerms );

            Complex qv(0, 0);
            Complex qw(0, 0);
            // Loop over the number of terms in the correlator
            for ( size_t i = 0; i < nTerms; i++ ){
                Complex tv;
                Complex tw;
                // loop over the # of products in each term (maximum of n)
                for ( size_t j = 0; j < lut[ i ].size(); j++ ){
                    NativeMask tm = maskAndCompactify( lut[i][j], m, 0, n );
                    
                    if ( 0 == tm ) continue;

                    std::bitset<MAX_SET_SIZE> bs( tm );
                    LOG_IF_F( INFO, DEBUG, "bs = %s", bs.to_string().c_str() );
                    
                    if ( qvm.count( bs ) == 0 ){
                        LOG_F( INFO, "NOT FOUND" );
                    }

                    auto q = qvm[ bs ];
                    LOG_IF_F( INFO, DEBUG, "t=%f + i%f", q.getQV().real(), q.getQV().imag() );

                    if ( 0 == j ){
                        tv = q.getQV();
                        tw = q.getW();
                    }
                    else {
                        tv *= q.getQV();
                        tw *= q.getW();
                    }
                }
                qv += tv;
                qw += tw;
                LOG_IF_F( INFO, DEBUG, "tv=%f + i%f", tv.real(), tv.imag() );
                LOG_IF_F( INFO, DEBUG, "tw=%f + i%f", tw.real(), tw.imag() );
            }
            this->v = qv;
            this->w = qw;

            LOG_F( INFO, "qv=%f + i%f", qv.real(), qv.imag() );
            LOG_F( INFO, "qw=%f + i%f", qw.real(), qw.imag() );
            LOG_F( INFO, "Finished computing n=%zu correlator", n );
        }

        Complex calculate(  ){
            return (this->v / this->w);
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

            LOG_IF_F( INFO, DEBUG, "%s = (im=%s, mm=%s, rm=%s)", std::bitset<8>( frm ).to_string().c_str(), std::bitset<8>(im).to_string().c_str(), std::bitset<8>(mm).to_string().c_str(), std::bitset<8>(rm).to_string().c_str() );
            
            return frm;
        } // maskAndCompactify
    };
}


#endif