#ifndef MULTICUMULANTS_TOYMCEVENT_H
#define MULTICUMULANTS_TOYMCEVENT_H

#include <string>
//ROOT include
#include <TH1.h>
//User
#include "ToyMC/ToyMCParticle.h"

/**
 * Class to define Sets and subsets
 * 
 */
namespace toymc
{

class ToyMCEvent
{
public:
	virtual const char* name() const { return "ToyMCEvent"; }
	virtual const char* classname() const { return "ToyMCEvent"; }

	/**
	 * Constructor - 
	 */
	ToyMCEvent()
        {
           set(500);
        }

	~ToyMCEvent()
        {
           std::vector<ToyMCParticle>().swap(this->_particles);
        }

        void set(int npart)
        {
           this->_particles.reserve(npart);
	}

        void add(ToyMCParticle part) 
        {
           this->_particles.push_back(part);
        };

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
        std::vector<ToyMCParticle> _particles;
};
}
#endif
