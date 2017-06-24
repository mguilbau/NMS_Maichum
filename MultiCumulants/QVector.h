#ifndef MULTICUMULANTS_QVECTOR_H
#define MULTICUMULANTS_QVECTOR_H

#include <iostream>
#include <string>
#include <bitset>
#include <unordered_map>
#include <vector>

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

		virtual void fill(double phi, double w )
		{
				int &power = this->_i;
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

		int _i;
		int _j;
protected:
		Complex _qvector;
		Complex _weights;
		Harmonic _harm;

};



typedef std::vector< std::vector<QVector> > QVectorVector;
// typedef std::vector< std::vector<bool> >    QVectorMask;


class QTerms {

public:
	void init_partition(std::vector<int> &s, std::vector<int> &m, int n, std::vector< std::bitset<MAX_SET_SIZE> > &bm, size_t &np ){

		s.resize( n );
		m.resize( n );
		bm.resize( n );

		for ( size_t i = 0; i < n; i++ ){
			s[i] = 1;
			m[i] = 1;
			bm[0].set( i );
		}

		np = 1;

	}

	int next_partition(std::vector<int> &s, std::vector<int> &m, int n, std::vector< std::bitset<MAX_SET_SIZE> > &bm, size_t &np ) {
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


		for ( size_t i = 0; i < n; i++ ) {
			bm[i].reset();
		}
		// converts [s] to bitmask
		for ( size_t i = 0; i < n; i++){
			bm[ s[i]-1 ].set( i );
		}

		np = max;

		return 1;
	}

	int factorial( int n ){
		if ( n <= 1 ) return 1;
		return n * factorial( n - 1 );
	}

	int coeff( int n ){
		return pow( -1, n-1 ) * factorial( n - 1 );
	}


	void test( size_t n ){
		LOG_F( INFO, "test(%d)", n );
		std::vector<int> s;
		std::vector<int> m;
		std::vector< std::bitset<MAX_SET_SIZE> > bm;
		size_t np = 0;
		init_partition( s, m, n, bm, np );

		size_t iTerm = 0;
		LOG_F( INFO, "np=%d [%d]", np, iTerm );
		LOG_F( INFO, "%d *\n%s", coeff(np), maskString( bm, np ).c_str() );


		while( next_partition( s, m, n, bm, np ) ){
			iTerm++;
			LOG_F( INFO, "np=%d [%d]", np, iTerm );
			LOG_F( INFO, "%d *\n%s", coeff(np), maskString( bm, np ).c_str() );
		}



		LOG_F( INFO, "test(%d) complete", n );
	}

	std::string maskString( std::vector< std::bitset<MAX_SET_SIZE> > bm, size_t np ){
		std::string s="";
		for ( size_t i = 0; i < np; i++){
			s+= "[" + bm[i].to_string() + "] \n";
		}
		return s;
	}

};



class QVectorMask {
public:
	std::bitset<MAX_SET_SIZE> bits;
	size_t i=0;
	size_t j=0;
};

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
			   this->generateBitmasks();

			   printf("# of Qvs=%d == %d\n", this->_qvm.size(), this->_masks.size() );

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


			// size_t nn = this->_masks.size();
			// for ( size_t i = 0; i < nn; i++ ){

			// 	if ( (this->_masks[i].bits & setMask) == this->_masks[i].bits ) {
			// 		// LOG_S(INFO) << "PASS : " << this->_masks[i].bits.to_string() << std::endl;
			// 		size_t &mi = this->_masks[i].i;
			// 		size_t &mj = this->_masks[i].j;
					
			// 		this->_q[mi][mj].fill(phi, weight, mi);
			// 	} else {
			// 		// LOG_S(INFO) << "FAIL : " << this->_masks[i].bits.to_string() << std::endl;
			// 	}
			// }
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
			std::string s="size=" + std::to_string( this->_masks.size() ) + "\n";
			for ( size_t i = 0; i < this->_masks.size(); i++){
				s+= "[" + std::to_string(this->_masks[i].i) + "]" + "[" + std::to_string(this->_masks[i].j) + "]\t" +  this->_masks[i].bits.to_string() + "\n";
			}
			return s;
		}

protected:
		Set _set;
		bool _useWeights;
		QVectorVector _q;
		std::vector<QVectorMask> _masks;
		std::unordered_map< std::bitset<MAX_SET_SIZE>, QVector > _qvm;
};



} // namespace cumulant
#endif
// Local Variables:
//  mode: C++
// End:
