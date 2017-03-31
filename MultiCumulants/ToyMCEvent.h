#ifndef MULTICUMULANTS_TOYMCEVENT_H
#define MULTICUMULANTS_TOYMCEVENT_H

#include <string>
//ROOT include
#include <TH1.h>
//User
#include "MultiCumulants/ToyMCParticle.h"

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
           set(500, 7, 0.04);
        }

	~ToyMCEvent()
        {
          std::vector<double>().swap(this->_vn);
        }

        void set(int npart, unsigned int nharm, double vn)
        {
           this->_particles.resize(npart);

           this->_vn.resize(nharm);
           for(unsigned int ih = 0; ih < this->_vn.size(); ++ih)
           {
             this->_vn[ih] = vn + 0.01*(ih+1);
           }
	}

        void add() {};

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
        std::vector<double> _vn;
};
}
#endif
