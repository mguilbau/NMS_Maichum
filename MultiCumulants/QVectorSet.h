#ifndef MULTICUMULANTS_QVECTOR_SET_H
#define MULTICUMULANTS_QVECTOR_SET_H


#include <iostream>
#include <string>
#include <bitset>
#include <unordered_map>
#include <vector>

#include "MultiCumulants/Types.h"
#include "MultiCumulants/Subsets.h"
#include "MultiCumulants/Algorithm.h"
#include "MultiCumulants/QVector.h"

// logging library
#define LOGURU_IMPLEMENTATION 1
#include "vendor/loguru/loguru.hpp"


namespace cumulant{

    typedef std::unordered_map< std::bitset<MAX_SET_SIZE>, QVector > QVectorMap;


    class QVectorSet{
    public:
        virtual const char* name() const { return "QVectorSet"; }
        virtual const char* classname() const { return "QVectorSet"; }

            //Constructors
            QVectorSet()
            : _set(0), _useWeights(false)
            {
            }

            QVectorSet(const HarmonicVector& h, Set set, bool useweights)
            {
                this->_set = set;
                if(this->_set.size() != h.size())  
                {
                    this->_set.resize(h.size());
                }

                this->_useWeights = useweights;
                this->generateBitmasks();

                LOG_F( INFO, "# of Qvs=%d == %d\n", this->_qvm.size(), this->_masks.size() );

            }

            //Destructors
            ~QVectorSet() {}

            virtual QVectorMap getQ() { return this->_qvm;}

            virtual void generateBitmasks()
            {
                algo::Combinations c;

                size_t n = this->_set.size();

                std::vector<int> ints;
                for (size_t i = 0; i < n; ints.push_back(i++));

                for(size_t k = 1; k <=n; ++k)
                {
                    size_t nC = 0;
                    do
                    {
                        std::bitset<MAX_SET_SIZE> bits;
                        QVectorMask mask;
                        mask.i=k-1;
                        mask.j=nC;
                        
                        for (size_t ik = 0; ik < k; ++ik)
                        {
                            mask.bits.set( ints[ik] );
                            bits.set( ints[ik] );
                        }
                        this->_masks.push_back( mask );
                        this->_qvm[ bits ]._i = k-1; // will create the kv pair and set the qv power
                        this->_qvm[ bits ]._j = nC;

                        ++nC;
                    } while(c.next_combination(ints.begin(), ints.begin() + k, ints.end()));
                }

            }

            virtual void fill(std::vector<double> &val, double &phi, double &w)
            {
                Real weight = (!this->_useWeights ? 1 : w);

                std::bitset<MAX_SET_SIZE> setMask = this->_set.setMask(val);

                for ( auto kv : this->_qvm ){
                    if ( ( kv.first & setMask ) == kv.first ){
                        kv.second.fill( phi, weight );
                    }
                }
            }

            virtual void reset()
            {
                for ( auto kv : this->_qvm ){
                    kv.second.reset(); 
                }
            }

            virtual std::string print()
            {
                std::string s = "";
                // for(size_t i = 0; i < this->_q.size(); ++i)
                // {
                //     for(size_t j = 0; j < this->_q[i].size(); ++j)
                //     {  
                //         s += "index (" + std::to_string(i) + ", " + std::to_string(j) + "): \n";
                //         s += this->_q[i][j].toString();
                //         s += "\n"; 
                //     }
                //     s += "\n";
                // }
                return s;
            }

            std::string toString(){
                    std::string s = "";

                    s += this->classname();
                    s += "< Set = " + this->_set.toString();
                    if(this->_useWeights) s += ", Use weights: TRUE ";
                    else                  s += ", Use weights: FALSE ";
                    s += ">";

                    return s;
            }

            std::string maskString(){
                std::string s="size=" + std::to_string( this->_masks.size() ) + "\n";
                for ( size_t i = 0; i < this->_masks.size(); i++){
                    s+= "[" + std::to_string(this->_masks[i].i) + "]" + "[" + std::to_string(this->_masks[i].j) + "]\t" +  this->_masks[i].bits.to_string() + "\n";
                }
                return s;
            }

    protected:
            Set _set;
            bool _useWeights;

            std::vector<QVectorMask> _masks;
            QVectorMap _qvm;
    };
} // napesapce cumulants


#endif