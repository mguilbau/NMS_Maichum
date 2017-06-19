#ifndef MULTICUMULANTS_CORRELATOR_H
#define MULTICUMULANTS_CORRELATOR_H

#include <iostream>
#include <stdio.h>
#include <vector>
#include <array>
#include <algorithm>
#include <math.h>
#include <set>
#include <string>
#include <cstring>
#include <map>
#include <stdlib.h>

#include "MultiCumulants/Types.h"
#include "MultiCumulants/QVector.h"

// logging library
#define LOGURU_IMPLEMENTATION 1
#include "vendor/loguru/loguru.hpp"

const int MAX_ORDER = 16;

namespace cumulant{
class Correlator
{
public:
	virtual const char* name() const { return "Correlator"; }
	virtual const char* classname() const {  return "Correlator"; }

	//Constructors
	Correlator()
	{
            LOG_S(INFO) << "Creating correlator container" << std::endl;
	}

	//Destructors
	~Correlator() {}

        //Function
        /*
            next
                - given the partitioning scheme represented by s and m, generate the next
            Returns: 1, if a valid partitioning was found 0, otherwise
        */
        int next(std::array<int, MAX_ORDER> &s, std::array<int, MAX_ORDER> &m, int n) {
            /* Update s: 1 1 1 1 -> 2 1 1 1 -> 1 2 1 1 -> 2 2 1 1 -> 3 2 1 1 -> 1 1 2 1 ... */
            int i = 0;
            ++s[i];
            while ((i < n - 1) && (s[i] > m[i] + 1)) {
                s[i] = 1;
                ++i;
                ++s[i];
            }
        
            /* If i has reached the n-1 element, then the last unique partitiong has been found*/
            if (i == n - 1)
                return 0;
        
            /* Because all the first i elements are now 1, s[i] (i + 1 element) is the largest. 
            So we update max by copying it to all the first i positions in m.*/
            int max = s[i];
            for (i = i - 1; i >= 0; --i)
                m[i] = max;
        
            return 1;
        }
        
        void print_as_Qvs( std::array<int, MAX_ORDER> s, int n, int np ){
        
            for ( int i = 1; i < np+1; i++ ){
                printf( "<" );
                char* space = "";
                for ( int j =0; j < n; j++ ){
                    if ( s[j] == i ){
                        printf( "%sQ_%d", space, j+1 );
                        space = " ";
                    }
                }
                printf( ">" );
            }
            
        }
        
        int max_in_map( std::map<int, int> &m ){
        
            int v = 0;
            for ( auto kv : m ){
                if ( kv.second > v ){
                    v = kv.second;
                }
            }
            return v;
        }
        
        int factorial( int n ){
            if ( n <= 1 ) return 1;
            return n * factorial( n - 1 );
        }
        int coeff( int n ){
            return pow( -1, n-1 ) * factorial( n - 1 );
        }

        void calculate(impl2::QVectorSet qvset) {}

protected:
};	

} // namespace cumulant
#endif
// Local Variables:
//  mode: C++
// End:
