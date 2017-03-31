#ifndef MULTICUMULANTS_TOYMCGENERATOR_H
#define MULTICUMULANTS_TOYMCGENERATOR_H

#include <string>

#include "MultiCumulants/ToyMCEvent.h"

/**
 * Class to define Sets and subsets
 * 
 */
namespace toymc
{

class ToyMCGenerator
{
public:
	virtual const char* name() const { return "ToyMCGenerator"; }
	virtual const char* classname() const { return "ToyMCGenerator"; }
	
	/**
	 * Constructor - 
	 */
	ToyMCGenerator(){}

	~ToyMCGenerator(){}

        void generate(){}
	
	/**
	 * builds the string representation of the object
	 *
	 * @return String representation of this
	 */
	std::string toString(){
		std::string s = "";
		s += this->classname(); 
		return s; 
	}

protected:
       
};
}
#endif
