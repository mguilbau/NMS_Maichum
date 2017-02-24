#ifndef MULTICUMULANTS_SUBSET_H
#define MULTICUMULANTS_SUBSET_H

#include "Types.h"

#include <string>



/**
 * Class to define an Eta Region
 * 
 */
class Subset
{
public:
	virtual const char* name() const { return "Subset"; }
	virtual const char* classname() const { return "Subset"; }
	
	/**
	 * Constructor - 
	 */
	Subset(){
                set( 0., 0., 0., 0., pCharge::kCharged );
        }

	Subset( Real etaMin, Real etaMax, Real ptMin, Real ptMax, pCharge charge ){
		set( etaMin, etaMax, ptMin, ptMax, charge );
	}

	~Subset(){}
	
        virtual void set( Real etaMin, Real etaMax, Real ptMin, Real ptMax, pCharge charge ){
                setEta(etaMin, etaMax);
                setPt(ptMin, ptMax);
                setCharge(charge);
        }

	virtual void setEta( Real etaMin, Real etaMax ){
		this->_etaMin = etaMin;
		this->_etaMax = etaMax;
	}

	virtual void setPt( Real ptMin, Real ptMax ){
		this->_ptMin = ptMin;
		this->_ptMax = ptMax;
	}

	virtual void setCharge( pCharge charge ){
		this->_charge = charge;
	}

	/**
	 * builds the string representation of the object
	 *
	 * @return String representation of this
	 */
	std::string toString(){
		std::string s = "";

		s += this->classname(); 
		s += " <min eta=" + std::to_string( this->_etaMin );
		s += ", max eta=" + std::to_string( this->_etaMax );
		s += ", min pt =" + std::to_string( this->_ptMin );
		s += ", max pt =" + std::to_string( this->_ptMax );
		s += ", charge =" + std::to_string( this->_charge );
		s += ">";

		return s; 
	}

protected:
	Real _etaMin;			/**< minimum eta **/ 
	Real _etaMax;			/**< maximum eta **/
	Real _ptMin;			/**< minimum pt **/ 
	Real _ptMax;			/**< maximum pt **/
        pCharge _charge;                /**< charge **/

};


#endif
