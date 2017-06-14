#ifndef MULTICUMULANTS_QVECTOR_H
#define MULTICUMULANTS_QVECTOR_H

#include <iostream>
#include <string>
#include <bitset>

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

namespace impl1{
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

			virtual void fill( std::vector<double> val, double phi, double w)
			{
				generateMask( val );
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

	protected:
			Set _set;
			bool _useWeights;
			QVectorVector _q;
			QVectorMask   _mask;
	};
} // namespace impl1


namespace impl2{

	typedef std::vector< std::vector<QVector> > QVectorVector;
	// typedef std::vector< std::vector<bool> >    QVectorMask;


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
						QVectorMask mask;
						mask.i=k-1;
						mask.j=nC;
						
						for (size_t ik = 0; ik < k; ++ik)
						{
							mask.bits.set( ints[ik] );
						}
						this->_masks.push_back( mask );
						++nC;
					} while(c.next_combination(ints.begin(), ints.begin() + k, ints.end()));
				}

			}

			virtual void fill(std::vector<double> &val, double &phi, double &w)
			{
				Real weight = (!this->_useWeights ? 1 : w);

				std::bitset<MAX_SET_SIZE> setMask = this->_set.setMask(val);
				// LOG_S(INFO) << "in Sets: " << setMask.to_string() << std::endl;

				size_t nn = this->_masks.size();
				for ( size_t i = 0; i < nn; i++ ){

					if ( (this->_masks[i].bits & setMask) == this->_masks[i].bits ) {
						// LOG_S(INFO) << "PASS : " << this->_masks[i].bits.to_string() << std::endl;
						size_t &mi = this->_masks[i].i;
						size_t &mj = this->_masks[i].j;
						
						this->_q[mi][mj].fill(phi, weight, mi);
					} else {
						// LOG_S(INFO) << "FAIL : " << this->_masks[i].bits.to_string() << std::endl;
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
				std::string s="size=" + std::to_string( this->_masks.size() ) + "\n";
				for ( size_t i = 0; i < this->_masks.size(); i++){
					s+= "[" + std::to_string(this->_masks[i].i) + "]" + "[" + std::to_string(this->_masks[i].j) + "]\t" +  this->_masks[i].bits.to_string() + "\n";
				}
				return s;
			}

			void make_coeff(){

				for ( size_t i = 0; i < this->_q.size(); i++ ){
					for(size_t j = 0; j < this->_q[i].size(); ++j){  
						
					}
				}


			}

	protected:
			Set _set;
			bool _useWeights;
			QVectorVector _q;
			std::vector<QVectorMask> _masks;
	};
} // namespace impl2

namespace impl3{

	typedef std::vector< std::vector<QVector> > QVectorVector;
	// typedef std::vector< std::vector<bool> >    QVectorMask;


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
						QVectorMask mask;
						mask.i=k-1;
						mask.j=nC;
						
						for (size_t ik = 0; ik < k; ++ik)
						{
							mask.bits.set( ints[ik] );
						}
						this->_masks.push_back( mask );
						++nC;
					} while(c.next_combination(ints.begin(), ints.begin() + k, ints.end()));
				}

			}

			virtual void fill(std::vector<double> &val, double &phi, double &w)
			{
				Real weight = (!this->_useWeights ? 1 : w);

				std::bitset<MAX_SET_SIZE> setMask = this->_set.setMask(val);
				// LOG_S(INFO) << "in Sets: " << setMask.to_string() << std::endl;

				size_t nn = this->_masks.size();
				if ( 15 == nn ){
					if ( (this->_masks[0].bits & setMask) == this->_masks[0].bits ) {
						size_t &mi = this->_masks[0].i;
						size_t &mj = this->_masks[0].j;
						
						this->_q[mi][mj].fill(phi, weight, mi);
					}
					if ( (this->_masks[1].bits & setMask) == this->_masks[1].bits ) {
						size_t &mi = this->_masks[1].i;
						size_t &mj = this->_masks[1].j;
						
						this->_q[mi][mj].fill(phi, weight, mi);
					}
					if ( (this->_masks[2].bits & setMask) == this->_masks[2].bits ) {
						size_t &mi = this->_masks[2].i;
						size_t &mj = this->_masks[2].j;
						
						this->_q[mi][mj].fill(phi, weight, mi);
					}
					if ( (this->_masks[3].bits & setMask) == this->_masks[3].bits ) {
						size_t &mi = this->_masks[3].i;
						size_t &mj = this->_masks[3].j;
						
						this->_q[mi][mj].fill(phi, weight, mi);
					}
					if ( (this->_masks[4].bits & setMask) == this->_masks[4].bits ) {
						size_t &mi = this->_masks[4].i;
						size_t &mj = this->_masks[4].j;
						
						this->_q[mi][mj].fill(phi, weight, mi);
					}
					if ( (this->_masks[5].bits & setMask) == this->_masks[5].bits ) {
						size_t &mi = this->_masks[5].i;
						size_t &mj = this->_masks[5].j;
						
						this->_q[mi][mj].fill(phi, weight, mi);
					}
					if ( (this->_masks[6].bits & setMask) == this->_masks[6].bits ) {
						size_t &mi = this->_masks[6].i;
						size_t &mj = this->_masks[6].j;
						
						this->_q[mi][mj].fill(phi, weight, mi);
					}
					if ( (this->_masks[7].bits & setMask) == this->_masks[7].bits ) {
						size_t &mi = this->_masks[7].i;
						size_t &mj = this->_masks[7].j;
						
						this->_q[mi][mj].fill(phi, weight, mi);
					}
					if ( (this->_masks[8].bits & setMask) == this->_masks[8].bits ) {
						size_t &mi = this->_masks[8].i;
						size_t &mj = this->_masks[8].j;
						
						this->_q[mi][mj].fill(phi, weight, mi);
					}
					if ( (this->_masks[9].bits & setMask) == this->_masks[9].bits ) {
						size_t &mi = this->_masks[9].i;
						size_t &mj = this->_masks[9].j;
						
						this->_q[mi][mj].fill(phi, weight, mi);
					}
					if ( (this->_masks[10].bits & setMask) == this->_masks[10].bits ) {
						size_t &mi = this->_masks[10].i;
						size_t &mj = this->_masks[10].j;
						
						this->_q[mi][mj].fill(phi, weight, mi);
					}
					if ( (this->_masks[11].bits & setMask) == this->_masks[11].bits ) {
						size_t &mi = this->_masks[11].i;
						size_t &mj = this->_masks[11].j;
						
						this->_q[mi][mj].fill(phi, weight, mi);
					}
					if ( (this->_masks[12].bits & setMask) == this->_masks[12].bits ) {
						size_t &mi = this->_masks[12].i;
						size_t &mj = this->_masks[12].j;
						
						this->_q[mi][mj].fill(phi, weight, mi);
					}
					if ( (this->_masks[13].bits & setMask) == this->_masks[13].bits ) {
						size_t &mi = this->_masks[13].i;
						size_t &mj = this->_masks[13].j;
						
						this->_q[mi][mj].fill(phi, weight, mi);
					}
					if ( (this->_masks[14].bits & setMask) == this->_masks[14].bits ) {
						size_t &mi = this->_masks[14].i;
						size_t &mj = this->_masks[14].j;
						
						this->_q[mi][mj].fill(phi, weight, mi);
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
				std::string s="";
				for ( size_t i = 0; i < this->_masks.size(); i++){
					s+= "[" + std::to_string(this->_masks[i].i) + "]" + "[" + std::to_string(this->_masks[i].j) + "] " +  this->_masks[i].bits.to_string() + "\n";
				}
				return s;
			}

	protected:
			Set _set;
			bool _useWeights;
			QVectorVector _q;
			std::vector<QVectorMask> _masks;
	};
} // namespace impl3


namespace impl4{

	typedef std::vector< std::vector<QVector> > QVectorVector;
	// typedef std::vector< std::vector<bool> >    QVectorMask;


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
						QVectorMask mask;
						mask.i=k-1;
						mask.j=nC;
						
						for (size_t ik = 0; ik < k; ++ik)
						{
							mask.bits.set( ints[ik] );
						}
						this->_masks.push_back( mask );
						++nC;
					} while(c.next_combination(ints.begin(), ints.begin() + k, ints.end()));
				}

			}

			virtual void fill(std::vector<double> &val, double &phi, double &w)
			{

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
				std::string s="";
				for ( size_t i = 0; i < this->_masks.size(); i++){
					s+= "[" + std::to_string(this->_masks[i].i) + "]" + "[" + std::to_string(this->_masks[i].j) + "] " +  this->_masks[i].bits.to_string() + "\n";
				}
				return s;
			}

	protected:
			Set _set;
			bool _useWeights;
			QVectorVector _q;
			std::vector<QVectorMask> _masks;
	};
} // namespace impl4

} // namespace cumulant
#endif
// Local Variables:
//  mode: C++
// End:
