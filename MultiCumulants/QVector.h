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

	class QVector {
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


	



	class QVectorMask {
	public:
		std::bitset<MAX_SET_SIZE> bits;
		size_t i=0;
		size_t j=0;
	};
} // namespace cumulant
#endif
// Local Variables:
//  mode: C++
// End:
