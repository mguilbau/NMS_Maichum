#ifndef MULTICUMULANTS_QVECTOR_H
#define MULTICUMULANTS_QVECTOR_H

#include <iostream>
#include <string>
#include "MultiCumulants/Types.h"
#include "MultiCumulants/Subsets.h"
#include "MultiCumulants/Algorithm.h"

// logging library
#define LOGURU_IMPLEMENTATION 1
#include "vendor/loguru/loguru.hpp"

namespace cumulant{
class QVector
{
public:
	virtual const char* name() const { return "QVector"; }
	virtual const char* classname() const {  return "QVector"; }

        //Constructors
        QVector()
          : _qvector(0), _weights(0), _harm(0)
        {
        }

        QVector(Harmonic harm)
          : _qvector(0), _weights(0), _harm(harm)
        {
        }

        //Destructors
        ~QVector() {}
	
        virtual const Complex& operator()(int index) const
        {
                if(index==0) return this->_weights;
                else         return this->_qvector;
        }

        virtual Complex& operator()(int index)
        {
                if(index==0) return this->_weights;
                else         return this->_qvector;
        }

        virtual QVector& operator*=(const QVector& qv)
        {
                this->_harm += qv._harm;
                this->_qvector *= qv(1);
                this->_weights *= qv(0);
                return *this;
        }

        virtual void setHarm(Harmonic harm) 
        {
                this->_harm = harm;
        }

        virtual Harmonic getHarm() 
        {
                return this->_harm;
        }

        virtual Complex getQV()
        {
               return this->_qvector;
        }

        virtual Complex getW()
        {
               return this->_weights;
        }

        virtual void fill(double phi, double w, int power)
        {
                this->_qvector += Complex(pow(w, power) * cos(static_cast<double>(this->_harm) * phi), 
                                          pow(w, power) * sin(static_cast<double>(this->_harm) * phi));
                this->_weights += Complex(pow(w, power) * cos(0.*phi), 
                                          pow(w, power) * sin(0.*phi));
        }

        virtual void reset(const Complex& q=Complex(0,0))
        {
                this->_qvector = q;
                this->_weights = q;
        }

        std::string toString(){
                std::string s = "";

                s += this->classname();
                s += "<harmonic =  " + std::to_string( this->_harm );
                s += "  qvector =  " + std::to_string( this->operator()(1).real() ) + " + " 
                                     + std::to_string( this->operator()(1).imag() ) + ".i";
                s += ", weights =  " + std::to_string( this->operator()(0).real() ) + " + " 
                                     + std::to_string( this->operator()(0).imag() ) + ".i";
                s += ">";

                return s;
        }
    
protected:
        Complex _qvector;
        Complex _weights;
        Harmonic _harm;	
};

typedef std::vector< std::vector<QVector> > QVectorVector;
typedef std::vector< std::vector<bool> >    QVectorMask;

class QVectorSet
{
public:
	virtual const char* name() const { return "QVectorSet"; }
	virtual const char* classname() const { return "QVectorSet"; }

        //Constructors
        QVectorSet()
         : _set(0), _useWeights(false), _q(0)
        {
        }

        QVectorSet(const HarmonicVector& h, Set set, bool useweights)
        {
               this->_set = set;
               if(this->_set.size() != h.size()) 
               {
                 this->_set.resize(h.size());
               }

               resize(h);

               this->_useWeights = useweights;
        }

        //Destructors
        ~QVectorSet() {}

        virtual QVectorVector getQ() { return this->_q;}

        virtual void resize(HarmonicVector h)
        {
                algo::Combinations c;

                size_t n = h.size();
                this->_q.resize(n);     

                std::vector<int> ints;
                for (size_t i = 0; i < n; ints.push_back(i++));

                for(size_t k = 1; k <=n; ++k)
                {
                  this->_q[k-1].resize(c.getCombinations(n,k), QVector());
                  size_t nC = 0;

                  do
                  {  
                     QVector qv;
                     for (size_t ik = 0; ik < k; ++ik)
                     {
                        qv *= QVector(h[ints[ik]]);
                     }
                     this->_q[k-1][nC] = qv;
                     ++nC;
                  }
                  while(c.next_combination(ints.begin(), ints.begin() + k, ints.end()));
                }
        }

        virtual void generateMask(std::vector<double> val)
        {
                algo::Combinations c;

                size_t n = this->_set.size();
                this->_mask.resize(n);     

                std::vector<int> ints;
                for (size_t i = 0; i < n; ints.push_back(i++));

                for(size_t k = 1; k <=n; ++k)
                {
                  this->_mask[k-1].resize(c.getCombinations(n,k), true);
                  size_t nC = 0;

                  do
                  {  
                     bool answer = true;
                     for (size_t ik = 0; ik < k; ++ik)
                     {
                        answer *= this->_set.isFromSet(val)[ints[ik]];
                     }
                     this->_mask[k-1][nC] = answer;
                     ++nC;
                  }
                  while(c.next_combination(ints.begin(), ints.begin() + k, ints.end()));
                }

        }

        virtual void fill(double phi, double w)
        {
               Real weight = (!this->_useWeights ? 1 : w);

               for(size_t i = 0; i < this->_q.size(); ++i)
               {
                 for(size_t j = 0; j < this->_q[i].size(); ++j)
                 {
                   if(this->_mask[i][j]) this->_q[i][j].fill(phi, weight, i); 
                 }
               }
        }

        virtual void reset()
        {
               for(size_t i = 0; i < this->_q.size(); ++i)
                 for(size_t j = 0; j < this->_q[i].size(); ++j)
                   this->_q[i][j].reset(); 
        }

        virtual std::string print()
        {
               std::string s = "";
               for(size_t i = 0; i < this->_q.size(); ++i)
               {
                 for(size_t j = 0; j < this->_q[i].size(); ++j)
                 {  
                    s += "index (" + std::to_string(i) + ", " + std::to_string(j) + "): \n";
                    s += this->_q[i][j].toString();
                    s += "\n"; 
                 }
                 s += "\n";
               }
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
            std::string s = "";
            size_t n = this->_mask.size();
            for ( size_t i = 0 ; i < n; i++ ){
                for ( size_t j = 0; j < this->_mask[i].size(); j++ ){
                    if ( this->_mask[i][j] )
                        s+= "[1]";
                    else 
                        s+= "[0]";
                }
                s+="\n";
            }
            return s;
        }

protected:
        Set _set;
        bool _useWeights;
        QVectorVector _q;
        QVectorMask   _mask;
};
}
#endif
// Local Variables:
//  mode: C++
// End:
