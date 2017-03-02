#ifndef MULTICUMULANTS_ETA_REGION_H
#define MULTICUMULANTS_ETA_REGION_H

#include "Types.h"

#include <string>

// logging library
#define LOGURU_IMPLEMENTATION 1
#include "loguru.hpp"

/**
 * Class to define an Eta Region
 * 
 */
class EtaRegion
{
public:
	virtual const char* name() const { return "Eta Region"; }
	virtual const char* classname() const { return "EtaRegion"; }
	
	/**
	 * Constructor - 
	 */
	EtaRegion() {}

	EtaRegion( Real min, Real max) {
		set( min, max );
	}
	~EtaRegion() {}
	

	virtual void set( Real min, Real max ){
		// You can log the scope of the current function call, this helps in profiling (not really needed here, just example)
		LOG_SCOPE_FUNCTION(INFO);
		this->_min = min;
		this->_max = max;
		LOG_F( INFO, "_min=%f, _max=%f", this->_min, this->_max );
	}

	/**
	 * builds the string representation of the object
	 *
	 * @return String representation of this
	 */
	std::string toString(){
		std::string s = "";

		s += this->classname(); 
		s += " <min=" + std::to_string( this->_min );
		s += ", max=" + std::to_string( this->_max );
		s += ">";

		return s; 
	}

protected:
	Real _min;			/**< minimum eta **/ 
	Real _max;			/**< maximum eta **/

};


#endif