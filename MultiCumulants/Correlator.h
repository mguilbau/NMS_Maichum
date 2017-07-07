#ifndef MULTICUMULANTS_CORRELATOR_H
#define MULTICUMULANTS_CORRELATOR_H

#include "MultiCumulants/QVectorSet.h"

#include "vendor/loguru/loguru.hpp"

namespace cumulant{

    class Correlator{

    public:
        Correlator() {

        }
        Correlator( size_t n, QVectorMap &qvm ){
            build( n, qvm );
        }

        void build( size_t n, QVectorMap &qvm ){
            LOG_F( INFO, "computing correlator for n=%zu", n );

            auto lut = NativeMaskLUTs[ n-3 ];

            
            size_t nTerms = lut.size();
            // size_t nTerms = nTermsLUT[ n-3 ];
            // // size_t nTerms = sizeof( QMasks[n-3] ) / sizeof( QMasks[n-3][0] );

            LOG_F( INFO, "nTerms = %zu", nTerms );

            Complex q(0, 0);
            // // Loop over the number of terms in the correlator
            for ( size_t i = 0; i < nTerms; i++ ){
                Complex t;
                // loop over the # of products in each term (maximum of n)
                for ( size_t j = 0; j < lut[ i ].size(); j++ ){
                    if ( 0 == lut[ i ][ j ] ) continue;
                    std::bitset<MAX_SET_SIZE> bs( lut[ i ][ j ] );
                    // LOG_F( INFO, "bs = %s", bs.to_string().c_str() );
                    if ( qvm.count( bs ) == 0 ){
                        LOG_F( INFO, "NOT FOUND" );
                    }

                    auto q = qvm[ bs ];
                    // LOG_F( INFO, "t=%f + i%f", q.getQV().real(), q.getQV().imag() );
                    if ( 0 == j )
                        t = q.getQV();
                    else 
                        t *= q.getQV();
                }
                q += t;
                // LOG_F( INFO, "t=%f + i%f", t.real(), t.imag() );
            }

            LOG_F( INFO, "q=%f + i%f", q.real(), q.imag() );
            LOG_F( INFO, "Finished computing n=%zu correlator", n );
        }
    };
}


#endif